
###############################################################################
# 
# Script title: Variable selection using least cost path iteration(s)
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
#' Reselect variables using LCP-refined covariates from a completed LCP
#' iteration instead of straight-line ones, to check whether the straight-
#' line-based selection in 02_select_variables.R still holds with LCPs
#'
#' Only used for the 4 EARLY subgroups -- selection never runs on
#' "late" subgroups
#'
#' Saves to results/fisher_variable_selection_lcp<iteration>.rds instead (e.g.
#' _lcp1.rds, _lcp2.rds -- each iteration's reselection run gets its own
#' file
#'
#' SELECTION_MODE -- "standard" or "robust": same toggle as
#' 02_select_variables.R. "standard" (default) runs select_variables_forced() at
#' one fixed hyperparameter setting. "robust" runs 
#' select_variables_robustness_check() (3 fixed settings) and keeps a candidate 
#' if >= robust_majority_threshold of them selected it.
#'
#' When COMPARE_TUNED_RMSE is TRUE, for any subgroup whose reselected
#' variable set differs from what's currently locked in, this fits a FULLY TUNED
#'  model using the reselected variables on the SAME LCP lines already used for
#' reselection above, and compares its pooled RMSE per lag against the
#' EXISTING tuned iteration-<iteration> model's RMSE (fit with whatever
#' variable set is currently locked in). 
#'
#' Even if the reselected variables win this comparison, adopting them would 
#' generate different LCP lines on the NEXT iteration going forward -- so a win 
#' here means "worth adopting and rerunning straight line and subsequent LCP
#' iterations to see if it holds up through its own LCP chain.
#' 
#' ALSO compares against the PREVIOUS iteration's (iteration - 1)
#' official model, not just the current iteration's -- RMSE doesn't
#' necessarily improve LCP-over-LCP (a later refinement can genuinely be
#' worse than the one before it).
#' ============================================================================

# Load required packages
sapply(c("doParallel", "purrr", "tibble", "dplyr"),
       require, character.only = T)

source("code/01_extract_covariates.R")
source("code/lcp_tools.R")       
source("code/selection_tools.R") 
source("code/cv_diagnostics.R")  

#' ---- User inputs ----
target_period <- "early"
lines_dir     <- "outputs/lcp_lines/"
model_dir     <- "models/"

response_col   <- "GD"
forced_vars    <- "euc_geog"
candidate_vars <- c("temp_min_early", "snowpack_early", "precip_early", "cc_lf", "dnbr", "burn_01", "burn_high_01", "tri_45_riley")

#' Can only pick variable from this group; if one is selected, the others are 
#' dropped as candidates from future models
exclusive_var_groups <- list(c("dnbr", "burn_01", "burn_high_01"))
lags_available <- c(1, 3, 5, 7, 9)   # the lags you actually ran iteration 1 for
top_n_lags     <- length(lags_available) - 2  # plain mean across all of them --
                                            
iteration      <- 1            # reselecting from this LCP iteration's lines
date_type      <- "max"
noCLS          <- 60           

SELECTION_MODE <- "robust"   # "standard" or "robust" 

robust_majority_threshold <- 2 # robust mode only: keep a candidate if
                                # selected by >= this many of the 3 fixed
                                # hyperparameter settings

COMPARE_TUNED_RMSE <- TRUE     
noCLS_tuned_fit <- 60          

#' ---------------------------------------------------------------

lines_path <- function(iteration, sex, region, lag, period, dir = lines_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_", lag, "_",
         date_type, "_", period, ".rds")
}

#' Matches 03_fit_straight_line.R/04_lcp_iteration_cycle.R's model_path()
#' 
model_path <- function(iteration, sex, region, lag, period, dir = model_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_", lag, "_",
         date_type, "_gbm_", period, ".rds")
}

lcp_reselect_cache_dir <- "results/lcp_reselection_extractions/"
dir.create(lcp_reselect_cache_dir, recursive = TRUE, showWarnings = FALSE)
lcp_reselect_cache_path <- function(label, iteration, dir = lcp_reselect_cache_dir) {
  paste0(dir, label, "_lcp", iteration, ".rds")
}

fold_results       <- readRDS("results/fisher_cv_folds.rds")
selection_straight <- readRDS(paste0("results/fisher_variable_selection_lcp1.rds"))
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>%
  dplyr::mutate(index = dplyr::row_number())

early_labels <- names(fold_results)[grepl(paste0("_", target_period, "$"), names(fold_results))]

#' AUTO-SKIP: drop any early label whose straight-line selection has no
#' landscape variable at all but only straight-line geographic distance
skipped_no_landscape_var <- early_labels %>%
  purrr::keep(~ !any(candidate_vars %in% selection_straight[[.x]]$selected_vars))

if (length(skipped_no_landscape_var) > 0) {
  message("Skipping (no landscape variable in straight-line selection, so no real LCP-", iteration,
          " lines exist to reselect from): ", paste(skipped_no_landscape_var, collapse = ", "))
}

early_labels <- setdiff(early_labels, skipped_no_landscape_var)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

selection_lcp1 <- early_labels %>%
  purrr::set_names() %>%
  purrr::map(function(label) {

    folds <- fold_results[[label]]
    if (is.null(folds)) return(NULL)

    sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
    tag       <- subgroup_tag(sub_pairs)
    period    <- target_period

    #' Catches a genuinely missing file for any reason (e.g. 04 hasn't been
    #' rerun yet for this subgroup) rather than erroring mid-run.
    lag_paths <- lags_available %>%
      purrr::set_names(as.character) %>%
      purrr::map_chr(~ lines_path(iteration, tag$sex, tag$region, .x, period))
    missing <- lag_paths[!file.exists(lag_paths)]
    if (length(missing) > 0) {
      message("Skipping ", label, " -- missing LCP-", iteration, " lines file(s): ",
              paste(missing, collapse = ", "))
      return(NULL)
    }

    time_varying_straight <- any(names(dynamic_var_settings) %in% selection_straight[[label]]$selected_vars)
    message("Reselecting from LCP-", iteration, ": ", label, " (", SELECTION_MODE, " mode)",
            if (!time_varying_straight)
              paste0(" -- NOTE: straight-line selection excluded all dynamic variables, so 04's shortcut means these lags' LCP-",
                     iteration, " lines are duplicates of one representative lag, not independent evidence")
            else "")

    sub_data_by_lag <- lags_available %>%
      purrr::set_names(as.character) %>%
      purrr::map(function(lag) {
        lcp_lines <- readRDS(lag_paths[[as.character(lag)]])

        sd <- extract_all_covariates(lcp_lines, sub_pairs, lag = lag, stacks = stacks,
                                      already_buffered = TRUE)
        stopifnot(nrow(sd) == folds$n_pairs)
        sd
      })

    saveRDS(sub_data_by_lag, lcp_reselect_cache_path(label, iteration))

    if (SELECTION_MODE == "robust") {

      robust <- select_variables_robustness_check(sub_data_by_lag, folds, forced_vars, candidate_vars,
                                                    response_col, top_n_lags = top_n_lags, exclusive_groups = exclusive_var_groups)

      message("  stability across ", length(default_gbm_param_settings), " hyperparameter settings:")
      print(robust$stability, n = Inf)

      keep <- robust$stability %>%
        dplyr::filter(n_settings_selected >= robust_majority_threshold) %>%
        dplyr::pull(variable)
      selected_vars <- c(forced_vars, candidate_vars[candidate_vars %in% keep])

      message("  majority-vote selection (>= ", robust_majority_threshold, "/",
              length(default_gbm_param_settings), " settings): ", paste(selected_vars, collapse = ", "))

      list(selected_vars = selected_vars, mode = "robust",
           history = robust$results_by_setting$modest$history,
           robust_results = robust$results_by_setting, stability = robust$stability)

    } else {
      res <- select_variables_forced(sub_data_by_lag, folds, forced_vars, candidate_vars, response_col,
                                      top_n_lags = top_n_lags, exclusive_groups = exclusive_var_groups)
      c(res, list(mode = "standard"))
    }
  })

parallel::stopCluster(cl)
gc()

dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(selection_lcp1, paste0("results/fisher_variable_selection_lcp", iteration, ".rds"))

#' ----------------------------------------------------------------------
#' Side-by-side comparison -- straight-line selection vs. LCP-1 reselection,
#' for whichever early subgroups actually got reselected above 
#' ----------------------------------------------------------------------
reselected_labels <- names(selection_lcp1) %>% purrr::keep(~ !is.null(selection_lcp1[[.x]]))

comparison <- reselected_labels %>%
  purrr::map_dfr(function(label) {
    res           <- selection_lcp1[[label]]
    straight_vars <- selection_straight[[label]]$selected_vars
    lcp_vars      <- res$selected_vars
    tibble::tibble(
      label          = label,
      mode           = if (is.null(res$mode)) "standard" else res$mode,
      straight_line  = paste(straight_vars, collapse = ", "),
      reselected     = paste(lcp_vars, collapse = ", "),
      changed        = !setequal(straight_vars, lcp_vars)
    )
  })


names(comparison)[names(comparison) == "reselected"] <- paste0("lcp", iteration)

print(comparison, n = Inf)

not_reselected <- setdiff(names(fold_results)[grepl(paste0("_", target_period, "$"), names(fold_results))], reselected_labels)
if (length(not_reselected) > 0) {
  message("Not reselected (see messages above for why): ", paste(not_reselected, collapse = ", "))
}

#' ----------------------------------------------------------------------
#' COMPARE_TUNED_RMSE -- Only for labels where changed == TRUE; fits the 
#' reselected variables with the SAME tuning procedure (fit_gbm_cv()) as the 
#' existing iteration-<iteration> model, on the SAME LCP lines, and compares 
#' pooled RMSE per lag.
#' ----------------------------------------------------------------------
tuned_rmse_comparison <- NULL

if (COMPARE_TUNED_RMSE) {

  labels_to_compare <- comparison %>% dplyr::filter(changed) %>% dplyr::pull(label)

  if (length(labels_to_compare) == 0) {
    message("\nCOMPARE_TUNED_RMSE: no subgroup's reselection changed the variable set -- nothing to compare.")
  } else {

    message("\nCOMPARE_TUNED_RMSE: fitting fully-tuned models for the reselected variable set(s) (",
            paste(labels_to_compare, collapse = ", "), ") against the existing tuned iteration-",
            iteration, " model(s)...")

    fit_cl <- parallel::makeCluster(noCLS_tuned_fit)
    doParallel::registerDoParallel(fit_cl)

    #' Reads whatever model is saved at a given iteration and
    #' returns its pooled RMSE shared by both the current-iteration ("old") and
    #' previous-iteration ("prev") lookups below
    official_rmse <- function(it, tag, lag, period, label) {
      path <- model_path(it, tag$sex, tag$region, lag, period)
      if (!file.exists(path)) {
        message("  ", label, " | lag ", lag, ": no existing iteration-", it,
                " model at ", path, " -- treating as NA.")
        return(NA_real_)
      }
      fit <- readRDS(path)
      compute_fold_metrics(dplyr::semi_join(fit$pred, fit$bestTune))$global$rmse
    }

    tuned_rmse_comparison <- tryCatch({
      labels_to_compare %>%
        purrr::map_dfr(function(label) {

          folds    <- fold_results[[label]]
          new_vars <- selection_lcp1[[label]]$selected_vars
          old_vars <- selection_straight[[label]]$selected_vars
          period   <- target_period
          tag      <- subgroup_tag(pairs %>% dplyr::filter(index %in% folds$orig_index))

          sub_data_by_lag <- readRDS(lcp_reselect_cache_path(label, iteration))

          lags_available %>%
            purrr::map_dfr(function(lag) {

              sub_data <- sub_data_by_lag[[as.character(lag)]]

              #' New fit: reselected variables, fully tuned, on these
              #' SAME LCP lines.
              new_fit  <- fit_gbm_cv(sub_data, new_vars, response_col, folds)
              new_rmse <- compute_fold_metrics(dplyr::semi_join(new_fit$pred, new_fit$bestTune))$global$rmse

              #' Existing model at this iteration -- whatever
              #' variable set is currently locked in (old_vars) -- AND
              #' the previous iteration's model, whatever
              #' variable set THAT happened to be fit with. RMSE isn't
              #' guaranteed to improve iteration-over-iteration (a later
              #' LCP refinement can genuinely be worse), so a reselected
              #' set only represents real forward progress if it beats
              #' both "old" and "previous".
              old_rmse  <- official_rmse(iteration,     tag, lag, period, label)
              prev_rmse <- official_rmse(iteration - 1, tag, lag, period, label)

              tibble::tibble(
                label = label, lag = lag,
                old_vars = paste(old_vars, collapse = ", "), new_vars = paste(new_vars, collapse = ", "),
                prev_iter = iteration - 1, prev_rmse = prev_rmse,
                old_iter = iteration, old_rmse = old_rmse,
                new_rmse = new_rmse,
                new_beats_old  = new_rmse < old_rmse,
                new_beats_prev = new_rmse < prev_rmse
              )
            })
        })
    }, finally = {
      parallel::stopCluster(fit_cl)
    })

    print(tuned_rmse_comparison, n = Inf)

    #' Mean of the best top_n_lags (not a plain mean across all of
    #' lags_available) 
    top_n_mean <- function(x, n) {
      x_sorted <- sort(x)   # sort() drops NA automatically
      mean(x_sorted[seq_len(min(n, length(x_sorted)))])
    }

    summary_by_label <- tuned_rmse_comparison %>%
      dplyr::group_by(label) %>%
      dplyr::summarise(n_lags_available = dplyr::n(),
                        prev_rmse_top_n = top_n_mean(prev_rmse, top_n_lags),
                        old_rmse_top_n  = top_n_mean(old_rmse, top_n_lags),
                        new_rmse_top_n  = top_n_mean(new_rmse, top_n_lags),
                        new_beats_old_top_n  = new_rmse_top_n < old_rmse_top_n,
                        new_beats_prev_top_n = new_rmse_top_n < prev_rmse_top_n,
                        worth_adopting = new_beats_old_top_n & new_beats_prev_top_n,
                        .groups = "drop")

    message("\nSummary (fully tuned, best ", top_n_lags, " of ", length(lags_available),
            " lags -- prev = iteration ", iteration - 1, "'s official model, old = iteration ",
            iteration, "'s official model (current locked-in vars), new = reselected vars on",
            " iteration ", iteration, "'s lines):")
    print(summary_by_label, n = Inf)
    message("worth_adopting requires beating BOTH prev and old -- beating old alone isn't enough",
            " if old is itself a regression from prev (RMSE doesn't always improve LCP-over-LCP).")
  }
}

as.data.frame(summary_by_label)
#' If changed == TRUE for a subgroup and you want to adopt the LCP-<iteration>
#' result going forward: copy that subgroup's entry from
#' selection_lcp<iteration> into results/fisher_variable_selection.rds (or just
#' overwrite the whole file if you're adopting it for all subgroups), then rerun
#' 03_fit_straight_line.R and 04_lcp_iteration_cycle.R for that subgroup(s)
#' so the new models reflect the updated variable set.

