
###############################################################################
# 
# Script title: Fitting straight-line GBM
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
## ============================================================================
## Fit iteration-0 (straight-line) GBMs: one per subgroup x lag.
## ============================================================================

# Load required packages
sapply(c("doParallel", "dplyr"),
       require, character.only = T)

source("code/01_extract_covariates.R")
source("code/fit_tools.R")

# ---- User inputs ----
response_col <- "GD"
date_type    <- "max"
lags         <- 1:10
outdir <- "models/"     # same folder LCP iterations save into --
                                # iteration 0 = straight-line, distinguished
                                # only by the "lcp_0_" filename prefix
noCLS        <- 60 # Number of cores for parallelization
# ---------------------------------------------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fold_results       <- readRDS("results/fisher_cv_folds.rds")
selection_results <- readRDS("results/fisher_variable_selection_lcp2.rds")
pairs <- readRDS("data/pairs_fisher_geomasked.rds") %>% 
  dplyr::mutate(index = dplyr::row_number())

stacks <- build_dynamic_stacks()

## ----------------------------------------------------------------------
extraction_cache_path <- function(label, dir = "results/straight_line_extractions/") paste0(dir, label, ".rds")

## If necessary, narrow to just the subgroup(s) whose variable set changed
subgroups_to_run <- names(fold_results)
## subgroups_to_run <- "female_south_early"

## ----------------------------------------------------------------

for (label in subgroups_to_run) {

  folds <- fold_results[[label]]
  sel   <- selection_results[[label]]
  if (is.null(folds) || is.null(sel)) {
    message("Skipping ", label, " (no folds/selection for this subgroup)")
    next
  }
  predictor_cols <- sel$selected_vars

  ## If no dynamic covariates are in this subgroup's selected set, this will
  ## fit the model for lag 1 and duplicate the fit for the other 9 lags
  time_varying <- any(names(dynamic_var_settings) %in% predictor_cols)

  message("Fitting straight-line: ", label, " (", length(lags), " lags",
          if (!time_varying) ", no time-varying predictor selected -- fitting once and duplicating across lags)"
          else ")")

  sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
  lines     <- lines_from_pairs(sub_pairs)

  ## NULL for the 4 late subgroups
  cache_path      <- extraction_cache_path(label)
  sub_data_by_lag <- if (file.exists(cache_path)) readRDS(cache_path) else NULL

  ## fresh cluster per subgroup 
  cl <- parallel::makeCluster(noCLS)
  doParallel::registerDoParallel(cl)

  fit_shared <- NULL   # only ever populated/reused when !time_varying

  for (lag in lags) {

    out_path <- paste0(outdir, "lcp_0_", tolower(unique(sub_pairs$sex)), "_",
                        unique(sub_pairs$region), "_lag_", lag, "_", date_type,
                        "_gbm_", sub("^.*_(early|late)$", "\\1", label), ".rds")

    if (!time_varying && !is.null(fit_shared)) {
      saveRDS(fit_shared, out_path)
      next
    }

    sub_data <- if (!is.null(sub_data_by_lag) && as.character(lag) %in% names(sub_data_by_lag)) {
      sub_data_by_lag[[as.character(lag)]]
    } else {
      extract_all_covariates(lines, sub_pairs, lag = lag, stacks = stacks)
    }

    ## sub_data must be in the same row order as folds$index/indexOut
    ## (positions into pairs_sub as built in 00_build_cv_folds.R)
    stopifnot(nrow(sub_data) == folds$n_pairs)

    fit <- fit_gbm_cv(sub_data, predictor_cols, response_col, folds)
    saveRDS(fit, out_path)

    if (!time_varying) fit_shared <- fit
  }

  parallel::stopCluster(cl)
  gc()   # encourage the master process to reclaim memory between subgroups too
}

