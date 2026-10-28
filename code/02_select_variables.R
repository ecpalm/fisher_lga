
###############################################################################
# 
# Script title: Variable selection
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

#============================================================================
#' Variable selection -- once per EARLY subgroup (e.g., female_north_early) on 
#' the iteration-0 (straight-line) data, evaluated across lags. 
#' The corresponding LATE subgroup inherits its early counterpart's variable set. 
#'
#' This code manually forces straight-line geographic distance to be retained
#' in every model and uses pooled/global RMSE (per lag) as the base criterion, 
#' averaged across user-chosen subset of lags. Stops adding variables
#' once none of the remaining candidates improve the across-lag average.
#'
#' GBM hyperparameters are FIXED (not tuned) during selection to avoid extremely
#' long runtimes.
#'
#' SELECTION_MODE -- "standard" or "robust":
#'   "standard": forward selection at ONE fixed hyperparameter setting 
#'   "robust": for each subgroup, runs select_variables_robustness_check()
#'     forward selection independently under 3 fixed hyperparameter
#'     settings (modest/generous/regularized; see default_gbm_param_settings
#'     in selection_tools.R). A candidate variable is retained in the final
#'     selected_vars if AT LEAST robust_majority_threshold of
#'     the 3 settings selected it.
#' ============================================================================

# Load required packages
sapply(c("doParallel", "purrr", "tibble", "dplyr"),
       require, character.only = T)

source("code/01_extract_covariates.R")
source("code/selection_tools.R")

#' ----------------------------------------------------------------------
#' Run for the 4 "early" subgroups only, across all 10 lags. 
#' ----------------------------------------------------------------------
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

#' ---- User inputs  ----
selection_path <- "results/fisher_variable_selection.rds"

response_col   <- "GD"
forced_vars    <- "euc_geog"
candidate_vars <- c("temp_min_early", "snowpack_early", "precip_early", "cc_lf", "dnbr", "burn_01", "burn_high_01", "tri_45_riley")

#' Can only pick variable from this group; if one is selected, the others are 
#' dropped as candidates from future models
exclusive_var_groups <- list(c("dnbr", "burn_01", "burn_high_01"))
lags           <- 1:10
top_n_lags     <- 3    # aggregate candidate set's RMSE over its best 3 lags
date_type      <- "max"
noCLS          <- 60   # change depending on your CPU's budget

SELECTION_MODE <- "robust"   # "standard" or "robust"

robust_majority_threshold <- 2 # robust mode only: keep a candidate if
                                # selected by >= this many of the 3 fixed
                                # hyperparameter settings
#' ---------------------------------------------------------------

fold_results <- readRDS("results/fisher_cv_folds.rds")
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>%
  dplyr::mutate(index = dplyr::row_number())

# Raster stacks for dynamic variables (annual layers)
stacks <- build_dynamic_stacks()

extraction_cache_dir <- "results/straight_line_extractions/"
dir.create(extraction_cache_dir, recursive = TRUE, showWarnings = FALSE)
extraction_cache_path <- function(label, dir = extraction_cache_dir) paste0(dir, label, ".rds")

subgroup_labels <- names(fold_results)
early_labels    <- subgroup_labels[grepl("_early$", subgroup_labels)]
late_labels     <- subgroup_labels[grepl("_late$", subgroup_labels)]

#' Define which early subgroups actually get run through selection 
subgroups_to_run <- early_labels

existing_selection_results <- if (file.exists(selection_path)) {
  readRDS(selection_path)
} else {
  NULL
}

#' Parallelization: one cluster for the entire run
cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

new_selection_results <- subgroups_to_run %>%
  purrr::set_names() %>%
  purrr::map(function(label) {

    message("Selecting variables: ", label, " (", SELECTION_MODE, " mode)")
    folds <- fold_results[[label]]
    if (is.null(folds)) return(NULL)

    sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
    lines     <- lines_from_pairs(sub_pairs)

    sub_data_by_lag <- lags %>%
      purrr::set_names(as.character) %>%
      purrr::map(function(lag) {
        sd <- extract_all_covariates(lines, sub_pairs, lag = lag, stacks = stacks)
        #' sub_data must stay row-aligned with folds$index/indexOut positions
        stopifnot(nrow(sd) == folds$n_pairs)
        sd
      })

    #' cache for 03_fit_straight_line.R 
    saveRDS(sub_data_by_lag, extraction_cache_path(label))

    if (SELECTION_MODE == "robust") {

      robust <- select_variables_robustness_check(sub_data_by_lag, folds, forced_vars, candidate_vars,
                                                    response_col, top_n_lags = top_n_lags, exclusive_groups = exclusive_var_groups)

      message("  stability across ", length(default_gbm_param_settings), " hyperparameter settings:")
      print(robust$stability, n = Inf)

      keep <- robust$stability %>%
        dplyr::filter(n_settings_selected >= robust_majority_threshold) %>%
        dplyr::pull(variable)
      #' preserve candidate_vars' original order 
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

#' Merge freshly-run results into whatever was already saved for other
#' labels (see subgroups_to_run above)
selection_results <- if (is.null(existing_selection_results)) list() else existing_selection_results
selection_results[names(new_selection_results)] <- new_selection_results

#' "late" subgroups inherit their early counterpart's result 
  for (late_label in late_labels) {
    early_label <- sub("_late$", "_early", late_label)
    if (early_label %in% names(new_selection_results)) {
      selection_results[[late_label]] <- selection_results[[early_label]]
    }
  }
selection_results <- selection_results[intersect(subgroup_labels, names(selection_results))]   # restore original label order

variable_selection <- names(selection_results) %>%
  purrr::map_dfr(function(label) {
    res <- selection_results[[label]]
    early_label <- sub("_late$", "_early", label)
    tibble::tibble(
      label = label,
      mode = if (is.null(res)) NA_character_ else if (is.null(res$mode)) "standard" else res$mode,
      selected_vars = if (is.null(res)) NA_character_ else paste(res$selected_vars, collapse = ", "),
      inherited_from = if (grepl("_late$", label)) early_label else NA_character_
    )
  })

print(variable_selection, n = Inf)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(selection_results, selection_path)

#' 03_fit_straight_line.R and 04_lcp_iteration_cycle.R both read this file
#' and pull selection_results[[label]]$selected_vars for that subgroup 
