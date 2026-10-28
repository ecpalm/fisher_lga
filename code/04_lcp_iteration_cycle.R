
###############################################################################
# 
# Script title: Least cost path cycle: fit GBM, create LCP, extract 
#               covariates, refit GBM
#
# Article title: Dynamic landscape genetics predicts disturbance-driven
#                declines in genetic connectivity in an endangered forest 
#                carnivore
#
# Article authors: Eric Palm, Erin Landguth, Katherine Zeller, Rebecca Green,
#                  Kristine Pilgrim, Sarah Stock, Craig Thompson, 
#                  Michael Schwartz, Jody Tucker
#
# Script author: Eric Palm
#
###############################################################################

#' ============================================================================
#' LCP iteration cycle -- for each time lag
#'   1. predict iteration (N-1)'s model onto each year's raster stack
#'   2. build LCP lines from those resistance surfaces 
#'   3. extract covariates along the new LCP lines 
#'   4. fit iteration N's model and save
#'
#'
#' NO-TIME-VARYING-PREDICTOR SHORTCUT: the dynamic variables are the only predictors that
#' vary with lag/year -- everything else build_prediction_stack() can
#' use (the 4 static rasters, the degenerate euc_geog constant layer) is
#' identical regardless of which year gets requested. So for a subgroup
#' whose selected variable set excludes all dynamic variables, EVERY year's resistance
#' surface at a given iteration is pixel-identical, which means
#' compute_lcp_lines() produces the same corridors regardless of lag, and
#' the resulting fit is the same regardless of lag too (same reasoning as
#' the shortcut added to 03_fit_straight_line.R, just propagated through
#' LCP line generation as well here). The driver below runs the real
#' n_iterations-round chain for exactly one representative lag in that case, then
#' duplicates each iteration's model AND lines file across the other 9
#' lags instead of recomputing LCP lines + refitting 9 more times per
#' iteration. Recomputed fresh from selection_results every run, so a
#' subgroup that regains a dynamic variable after a future reselection just goes back
#' to running all 10 lags automatically.
#'
#' RESUMABILITY: before computing iteration N for a (subgroup, lag), the
#' driver checks whether model_path(N, ...) and lines_path(N, ...)
#' already exist and skips straight to iteration N+1 if so (unless
#' FORCE_RECOMPUTE is TRUE). 
#'
#' SUBGROUPS_TO_RUN: defaults to every subgroup (unchanged behavior).
#' Narrow this to just the subgroup(s) whose variable set actually
#' changed to leave every other subgroup's models/lcp_lines files
#' completely untouched -- combine with FORCE_RECOMPUTE <- TRUE for
#' exactly that workflow: subgroups_to_run scopes WHICH labels this run
#' even looks at (everything else is skipped before RESUMABILITY ever
#' gets a chance to apply), and FORCE_RECOMPUTE ensures the labels that
#' ARE in scope get a genuine recompute rather than resumability skipping
#' their now-stale iteration 1+ files. REMEMBER: an early subgroup's
#' "late" counterpart inherited the same reselected set verbatim in
#' 02_select_variables.R -- both belong in subgroups_to_run, not just
#' the early one.
#' ============================================================================
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# Load required packages
sapply(c("sf", "caret", "doParallel", "purrr", "dplyr"),
       require, character.only = T)

source("code/lcp_tools.R")

#' ---- User config -- adjust to your actual paths/columns ----
response_col  <- "GD"
date_type     <- "max"
model_dir     <- "models/"      # ALL models -- straight-line (iteration 0)
                                 # and LCP (iterations 1-3 by default) alike,
                                 # one folder
lines_dir     <- "outputs/lcp_lines/"
n_iterations  <- 1            # number of LCP refinement rounds after the
                                 # straight-line model 
FORCE_RECOMPUTE <- TRUE        # if TRUE, recomputes every iteration from scratch
                                 # even if its output files already exist 
lags          <- 1:10
representative_lag <- lags[1]   # used for models with only static predictors
noCLS_extract <- 100            # LCP line computation + covariate extraction
noCLS_model   <- 60             # caret::train() grid search -- heavier per-worker
#' ---------------------------------------------------------------

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lines_dir, recursive = TRUE, showWarnings = FALSE)

fold_results      <- readRDS("results/fisher_cv_folds.rds")
selection_results <- readRDS(paste0("results/fisher_variable_selection_final.rds"))
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>% 
  dplyr::mutate(index = dplyr::row_number())

#' ----------------------------------------------------------------------
#' iteration 0 = straight-line, 1-3 = LCP refinements -- all saved into the
#' same model_dir now, distinguished only by the "lcp_{iteration}_" prefix.
#' lines_path() mirrors it for the LCP line files run_lcp_iteration() saves.
#' ----------------------------------------------------------------------
model_path <- function(iteration, sex, region, lag, period, dir = model_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_", lag, "_",
         date_type, "_gbm_", period, ".rds")
}

lines_path <- function(iteration, sex, region, lag, period, dir = lines_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_", lag, "_",
         date_type, "_", period, ".rds")
}

#' ----------------------------------------------------------------------
#' One full iteration for one subgroup x lag. Manages its OWN two
#' clusters internally (unlike compute_lcp_lines()/extract_all_
#' covariates(), which never do) -- see the TWO CLUSTER SIZES note near
#' the top of this file for why extract_cl and model_cl each get created,
#' used, and fully torn down in sequence rather than the caller handing
#' in pre-made cluster objects.
#' ----------------------------------------------------------------------
run_lcp_iteration <- function(label, lag, iteration) {

  message("LCP iteration ", iteration, ": ", label, " | lag ", lag)

  folds <- fold_results[[label]]
  sel   <- selection_results[[label]]
  if (is.null(folds) || is.null(sel)) return(invisible(NULL))
  predictor_cols <- sel$selected_vars

  period <- sub("^.*_(early|late)$", "\\1", label)
  yr <- period_years %>% dplyr::filter(period == .env$period)

  sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
  domain    <- folds$model_domain
  tag       <- subgroup_tag(sub_pairs)

  prev_fit <- readRDS(model_path(iteration - 1, tag$sex, tag$region, lag, period))

  euc_geog_value <- stats::median(sub_pairs$euc_geog)
  years <- (yr$year_min - lag):(yr$year_max - lag)

  resistance_by_year <- years %>%
    purrr::set_names(as.character) %>%
    purrr::map(~ predict_resistance_surface(prev_fit, predictor_cols, .x, euc_geog_value, domain))

  #' Extraction phase: own cluster, guaranteed torn down (tryCatch's
  #' finally =) before the function moves on -- see the TWO CLUSTER SIZES
  #' note near the top of this file for why this can't just be a
  #' registerDoParallel() switch against an already-alive model_cl.
  extract_cl <- parallel::makeCluster(noCLS_extract)
  parallel::clusterEvalQ(extract_cl, {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    NULL
  })
  doParallel::registerDoParallel(extract_cl)

  sub_data <- tryCatch({
    lcp_lines <- compute_lcp_lines(resistance_by_year, sub_pairs, lag, domain)
    saveRDS(lcp_lines, lines_path(iteration, tag$sex, tag$region, lag, period))

    #' already_buffered = TRUE -- lcp_lines came out of compute_lcp_lines()
    #' above already buffered 500m and water-erased; extract_all_covariates()
    #' must NOT buffer it again (see the note on extract_static_covariates()
    #' in 01_extract_covariates.R).
    sd <- extract_all_covariates(lcp_lines, sub_pairs, lag = lag, stacks = stacks,
                                  already_buffered = TRUE)
    stopifnot(nrow(sd) == folds$n_pairs)
    sd
  }, finally = {
    parallel::stopCluster(extract_cl)   # fully closed before model_cl exists
  })

  #' Model-fitting phase: separate cluster, only created now that
  #' extract_cl's connections are freed.
  model_cl <- parallel::makeCluster(noCLS_model)
  parallel::clusterEvalQ(model_cl, {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    NULL
  })
  doParallel::registerDoParallel(model_cl)

  fit <- tryCatch({
    fit_gbm_cv(sub_data, predictor_cols, response_col, folds)
  }, finally = {
    parallel::stopCluster(model_cl)
  })

  saveRDS(fit, model_path(iteration, tag$sex, tag$region, lag, period))
}

#' ----------------------------------------------------------------------
#' Run all 8 subgroups x 10 lags, iterations 1..n_iterations IN ORDER PER
#' (subgroup, lag) -- iteration N only ever reads the PREVIOUS iteration's
#' model for that SAME (label, lag) (see model_path(iteration - 1, ...)
#' inside run_lcp_iteration() above)
#' ----------------------------------------------------------------------
all_failures <- character(0)

#' ---- ADAPT: narrow to just the subgroup(s) whose variable set changed
#' -- see the SUBGROUPS_TO_RUN header note above. Pair with
#' FORCE_RECOMPUTE <- TRUE when doing this after adopting a reselected
#' variable set. ----
subgroups_to_run <- names(fold_results)
#subgroups_to_run <- "female_south_early"
#' ----------------------------------------------------------------

for (label in subgroups_to_run) {

  folds <- fold_results[[label]]
  sel   <- selection_results[[label]]
  if (is.null(folds) || is.null(sel)) {
    message("Skipping ", label, " (no folds/selection for this subgroup)")
    next
  }
  predictor_cols    <- sel$selected_vars
  time_varying      <- any(names(dynamic_var_settings) %in% predictor_cols)
  has_landscape_var <- any(c(names(static_rasters), names(dynamic_var_settings)) %in% predictor_cols)

  sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
  tag       <- subgroup_tag(sub_pairs)
  period    <- sub("^.*_(early|late)$", "\\1", label)

  #' ------------------------------------------------------------------
  #' NO LANDSCAPE VARIABLE AT ALL -- see the NO-LANDSCAPE-VARIABLE-AT-ALL
  #' SHORTCUT note near the top of this file. Never calls
  #' run_lcp_iteration() for this subgroup; duplicates iteration 0's
  #' model and the straight-line geometry across every iteration and lag
  #' instead, then moves on to the next subgroup.
  #' ------------------------------------------------------------------
  if (!has_landscape_var) {
    message("LCP iterations: ", label, " -- no landscape variable selected (predictor_cols = ",
            paste(predictor_cols, collapse = ", "), "); the resistance surface would be spatially",
            " uniform, so LCP refinement is a no-op here. Duplicating the straight-line model",
            " across ", n_iterations, " iteration(s) and all ", length(lags), " lags instead of",
            " attempting LCP line computation.")

    lines0 <- lines_from_pairs(sub_pairs) %>% sf::st_buffer(500)

    for (lag in lags) {
      fit0_path <- model_path(0, tag$sex, tag$region, lag, period)
      if (!file.exists(fit0_path)) {
        message("  ", label, " | lag ", lag, ": no straight-line model at ", fit0_path,
                " -- skipping (run 03_fit_straight_line.R for this lag first).")
        all_failures <- c(all_failures, paste0(label, " | lag ", lag, " (missing iteration-0 model)"))
        next
      }
      fit0 <- readRDS(fit0_path)
      for (iteration in seq_len(n_iterations)) {
        saveRDS(fit0,   model_path(iteration, tag$sex, tag$region, lag, period))
        saveRDS(lines0, lines_path(iteration, tag$sex, tag$region, lag, period))
      }
    }

    next   # skip the time_varying/run_lcp_iteration logic below entirely for this subgroup
  }

  lags_to_run <- if (time_varying) lags else representative_lag

  message("LCP iterations: ", label,
          if (!time_varying)
            paste0(" -- no time-varying predictor selected, running only lag ", representative_lag,
                   " and duplicating across the other ", length(lags) - 1, " lags")
          else "")

  failed_lags <- character(0)

  for (lag in lags_to_run) {
    tryCatch({
      for (iteration in seq_len(n_iterations)) {

        already_done <- !FORCE_RECOMPUTE &&
          file.exists(model_path(iteration, tag$sex, tag$region, lag, period)) &&
          file.exists(lines_path(iteration, tag$sex, tag$region, lag, period))

        if (already_done) {
          message("  ", label, " | lag ", lag, " | iteration ", iteration,
                  ": already computed, skipping (set FORCE_RECOMPUTE <- TRUE to redo).")
          next
        }

        run_lcp_iteration(label, lag, iteration)
      }
    }, error = function(e) {
      message("FAILED: ", label, " | lag ", lag, " -- ", conditionMessage(e))
      message("  Skipping this (subgroup, lag); continuing with the rest of the run.")
      failed_lags <<- c(failed_lags, as.character(lag))
      all_failures <<- c(all_failures, paste0(label, " | lag ", lag))
    })
    gc()   # encourage the master process to reclaim memory between (subgroup, lag) combos too
  }

  if (!time_varying) {
    if (as.character(representative_lag) %in% failed_lags) {
      message("  ", label, ": representative_lag ", representative_lag,
              " failed -- nothing to duplicate across the other lags, skipping that step for this subgroup.")
    } else {
      #' representative_lag's iteration chain succeeded; duplicate its
      #' model + lines files across the other lags per iteration --
      #' exact, not approximate, per the header note above.
      other_lags <- setdiff(lags, representative_lag)
      for (iteration in seq_len(n_iterations)) {
        fit_shared   <- readRDS(model_path(iteration, tag$sex, tag$region, representative_lag, period))
        lines_shared <- readRDS(lines_path(iteration, tag$sex, tag$region, representative_lag, period))
        for (lag in other_lags) {
          saveRDS(fit_shared,   model_path(iteration, tag$sex, tag$region, lag, period))
          saveRDS(lines_shared, lines_path(iteration, tag$sex, tag$region, lag, period))
        }
      }
    }
  }
}

if (length(all_failures) > 0) {
  message("\n=== ", length(all_failures), " (subgroup, lag) combination(s) failed this run ===")
  message(paste("  -", all_failures, collapse = "\n"))
  message("Check models/ and outputs/lcp_lines/ for these before treating the run as complete.")
} else {
  message("\nAll (subgroup, lag) combinations completed without error.")
}

