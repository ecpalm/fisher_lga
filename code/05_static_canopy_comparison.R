
###############################################################################
# 
# Script title: Least cost path cycle with static vegetation variables for comparison 
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
#' Static vegetation comparison -- straight-line + 5 LCP iterations per
#' subgroup, reusing each subgroup's selected variable set but replacing the 
#' dynamic, lag-varying extraction of dynamic variables with a single fixed-year 
#' snapshot: 2011 for Early subgroups, 2019 for Late subgroups. 
#' ============================================================================

# Load required packages
sapply(c("doParallel", "tibble", "dplyr"),
       require, character.only = T)

source("code/lcp_tools.R")

#' ---- User configure ----
response_col    <- "GD"
date_type       <- "max"
model_dir       <- "models/"
lines_dir       <- "outputs/lcp_lines/"
n_iterations    <- 5
FORCE_RECOMPUTE <- FALSE
noCLS_extract   <- 100
noCLS_model     <- 60

static_years <- tibble::tibble(period = c("early", "late"), year = c(2011, 2019))

# ---------------------------------------------------------------

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lines_dir, recursive = TRUE, showWarnings = FALSE)

fold_results      <- readRDS("results/fisher_cv_folds.rds")

# Final variable selection used for top LCP models
selection_results <- readRDS("results/fisher_variable_selection_lcp2.rds") 
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>% 
  dplyr::mutate(index = dplyr::row_number())

model_path_static <- function(iteration, sex, region, year, period, dir = model_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_static_", year, "_",
         date_type, "_gbm_", period, ".rds")
}
lines_path_static <- function(iteration, sex, region, year, period, dir = lines_dir) {
  paste0(dir, "lcp_", iteration, "_", sex, "_", region, "_lag_static_", year, "_",
         date_type, "_", period, ".rds")
}

all_failures <- character(0)

for (label in names(fold_results)) {

  folds <- fold_results[[label]]
  sel   <- selection_results[[label]]
  if (is.null(folds) || is.null(sel)) {
    message("Skipping ", label, " (no folds/selection for this subgroup)")
    next
  }
  predictor_cols <- sel$selected_vars

  if (!any(names(dynamic_var_settings) %in% predictor_cols)) {
    message("Skipping ", label, " -- no dynamic variable in the selected variable set (",
            paste(predictor_cols, collapse = ", "),
            "); a static run would be identical to the existing dynamic run. Nothing to compare.")
    next
  }

  period <- sub("^.*_(early|late)$", "\\1", label)
  year   <- static_years$year[static_years$period == period]

  sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
  tag       <- subgroup_tag(sub_pairs)
  domain    <- folds$model_domain
  has_landscape_var <- any(c(names(static_rasters), names(dynamic_var_settings)) %in% predictor_cols)

  message("Static vegetation comparison: ", label, " (year ", year, ")")

  tryCatch({

    #' ---- Iteration 0: straight-line, static year ----
    fit0_path <- model_path_static(0, tag$sex, tag$region, year, period)

    if (FORCE_RECOMPUTE || !file.exists(fit0_path)) {
      lines0 <- lines_from_pairs(sub_pairs)

      extract_cl <- parallel::makeCluster(noCLS_extract)
      doParallel::registerDoParallel(extract_cl)
      sub_data0 <- tryCatch(
        extract_all_covariates_static(lines0, sub_pairs, year = year, stacks = stacks, vars = names(stacks)),
        finally = parallel::stopCluster(extract_cl)
      )
      stopifnot(nrow(sub_data0) == folds$n_pairs)

      model_cl <- parallel::makeCluster(noCLS_model)
      doParallel::registerDoParallel(model_cl)
      fit0 <- tryCatch(
        fit_gbm_cv(sub_data0, predictor_cols, response_col, folds),
        finally = parallel::stopCluster(model_cl)
      )
      saveRDS(fit0, fit0_path)
    } else {
      message("  iteration 0: already computed, skipping (set FORCE_RECOMPUTE <- TRUE to redo).")
    }

      #' ---- Iterations 1-5: LCP refinement using the fixed-year surface ----
      euc_geog_value <- stats::median(sub_pairs$euc_geog)

      for (iteration in seq_len(n_iterations)) {

        already_done <- !FORCE_RECOMPUTE &&
          file.exists(model_path_static(iteration, tag$sex, tag$region, year, period)) &&
          file.exists(lines_path_static(iteration, tag$sex, tag$region, year, period))

        if (already_done) {
          message("  iteration ", iteration, ": already computed, skipping",
                  " (set FORCE_RECOMPUTE <- TRUE to redo).")
          next
        }

        message("  static LCP iteration ", iteration, " (year ", year, ")")

        prev_fit <- readRDS(model_path_static(iteration - 1, tag$sex, tag$region, year, period))
        resistance <- predict_resistance_surface(prev_fit, predictor_cols, year, euc_geog_value, domain)

        extract_cl <- parallel::makeCluster(noCLS_extract)
        doParallel::registerDoParallel(extract_cl)
        sub_data <- tryCatch({
          lcp_lines <- compute_lcp_lines_static(resistance, sub_pairs, domain)
          saveRDS(lcp_lines, lines_path_static(iteration, tag$sex, tag$region, year, period))

          sd <- extract_all_covariates_static(lcp_lines, sub_pairs, year = year, stacks = stacks, vars = names(stacks),
                                               already_buffered = TRUE)
          stopifnot(nrow(sd) == folds$n_pairs)
          sd
        }, finally = parallel::stopCluster(extract_cl))

        model_cl <- parallel::makeCluster(noCLS_model)
        doParallel::registerDoParallel(model_cl)
        fit <- tryCatch(
          fit_gbm_cv(sub_data, predictor_cols, response_col, folds),
          finally = parallel::stopCluster(model_cl)
        )
        saveRDS(fit, model_path_static(iteration, tag$sex, tag$region, year, period))
        
      }

  }, error = function(e) {
    message("FAILED: ", label, " -- ", conditionMessage(e))
    message("  Skipping this subgroup; continuing with the rest of the run.")
    all_failures <<- c(all_failures, label)
  })

  gc()
}

if (length(all_failures) > 0) {
  message("\n=== ", length(all_failures), " subgroup(s) failed this run ===")
  message(paste("  -", all_failures, collapse = "\n"))
} else {
  message("\nAll eligible subgroups completed without error.")
}
