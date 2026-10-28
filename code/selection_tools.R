
###############################################################################
# 
# Script title: GBM variable selection tools
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
## Forward-selection tools 
## ============================================================================

# Load required packages
sapply(c("gbm", "CAST", "purrr", "tibble", "dplyr"),
       require, character.only = T)

## ----------------------------------------------------------------------
## Fit and evaluate one candidate variable set against one subgroup's
## folds, for ONE lag. Returns pooled RMSE across all held-out rows across
## all folds as in CAST::global_validation() 
##
## ----------------------------------------------------------------------
## gbm::gbm.fit() won't fit a model if there are too few observations 
## within a training fold. Caps it ONCE per subgroup from
## the smallest training fold so every fold uses the same hyperparameters 
cap_minobsinnode <- function(n_minobsinnode, folds, bag_fraction = 0.5) {
  min_train <- min(lengths(folds$index))
  cap <- max(1, floor((bag_fraction * min_train - 1) / 2 - 1e-9))
  min(n_minobsinnode, cap)
}

evaluate_variable_set <- function(sub_data, folds, vars, response_col,
                                   gbm_params = list(n.trees = 1000, interaction.depth = 1,
                                                      shrinkage = 0.01, n.minobsinnode = 10),
                                   seed = 1234) {

  set.seed(seed)
  gbm_params$n.minobsinnode <- cap_minobsinnode(gbm_params$n.minobsinnode, folds)

  fold_formula <- stats::reformulate(sprintf("`%s`", vars), response = sprintf("`%s`", response_col))

  fold_preds <- names(folds$index) %>%
    purrr::map_dfr(function(fold_name) {
      train_data <- sub_data[folds$index[[fold_name]], c(vars, response_col)]
      test_data  <- sub_data[folds$indexOut[[fold_name]], c(vars, response_col)]

      fold_fit <- gbm::gbm(
        formula = fold_formula, data = train_data, distribution = "gaussian",
        n.trees = gbm_params$n.trees, interaction.depth = gbm_params$interaction.depth,
        shrinkage = gbm_params$shrinkage, n.minobsinnode = gbm_params$n.minobsinnode
      )

      pred <- gbm::predict.gbm(fold_fit, newdata = test_data, n.trees = gbm_params$n.trees)
      tibble::tibble(obs = test_data[[response_col]], pred = pred)
    })

  sqrt(mean((fold_preds$obs - fold_preds$pred)^2))
}

## Same, but aggregated across lags as the mean of the BEST top_n_lags (of
## however many are in sub_data_by_lag) pooled RMSE values, not a plain
## mean across all of them. 
##
## top_n_lags is picked independently per candidate set (not a single
## fixed top-N list applied to all candidates) 
evaluate_variable_set_multilag <- function(sub_data_by_lag, folds, vars, response_col,
                                            gbm_params = list(n.trees = 1000, interaction.depth = 1,
                                                               shrinkage = 0.01, n.minobsinnode = 10),
                                            top_n_lags = 5, seed = 1234) {
  rmse_by_lag <- sub_data_by_lag %>%
    purrr::map_dbl(~ evaluate_variable_set(.x, folds, vars, response_col, gbm_params, seed))
  mean(sort(rmse_by_lag)[seq_len(min(top_n_lags, length(rmse_by_lag)))])
}

## ----------------------------------------------------------------------
## Forward selection with a mandatory base set (euc_geog), criterion =
## mean pooled RMSE across each candidate's best top_n_lags lags. 
## ----------------------------------------------------------------------
select_variables_forced <- function(sub_data_by_lag, folds, forced_vars, candidate_vars, response_col,
                                     gbm_params = list(n.trees = 1000, interaction.depth = 1,
                                                        shrinkage = 0.01, n.minobsinnode = 10),
                                     top_n_lags = 5, seed = 1234, exclusive_groups = list()) {

  eff_mino <- cap_minobsinnode(gbm_params$n.minobsinnode, folds)
  if (eff_mino < gbm_params$n.minobsinnode) {
    message("  NOTE: n.minobsinnode capped ", gbm_params$n.minobsinnode, " -> ", eff_mino,
            " (smallest training fold has ", min(lengths(folds$index)),
            " rows; gbm needs nTrain * bag.fraction > 2 * n.minobsinnode + 1)")
  }

  drop_groupmates <- function(added_var, remaining) {
    hit <- exclusive_groups %>% purrr::keep(~ added_var %in% .x)
    if (length(hit) == 0) return(remaining)
    excluded <- intersect(setdiff(unlist(hit), added_var), remaining)
    if (length(excluded) > 0) {
      message("  ", added_var, " is in a mutually-exclusive group -- removing ",
              paste(excluded, collapse = ", "), " from further consideration")
    }
    setdiff(remaining, excluded)
  }

  current_vars <- forced_vars
  remaining    <- candidate_vars
  for (fv in forced_vars) remaining <- drop_groupmates(fv, remaining)
  current_rmse <- evaluate_variable_set_multilag(sub_data_by_lag, folds, current_vars, response_col,
                                                  gbm_params, top_n_lags, seed)

  history <- tibble::tibble(step = 0, added = NA_character_,
                             vars = paste(current_vars, collapse = ", "), mean_rmse_top_lags = current_rmse)

  step <- 1
  repeat {
    if (length(remaining) == 0) break

    trial_rmse <- remaining %>%
      purrr::map_dbl(~ evaluate_variable_set_multilag(sub_data_by_lag, folds, c(current_vars, .x),
                                                        response_col, gbm_params, top_n_lags, seed))

    best_i <- which.min(trial_rmse)

    if (trial_rmse[best_i] >= current_rmse) break   # early stopping: no remaining candidate improves RMSE

    current_vars <- c(current_vars, remaining[best_i])
    current_rmse <- trial_rmse[best_i]
    history <- history %>%
      dplyr::bind_rows(tibble::tibble(step = step, added = remaining[best_i],
                                       vars = paste(current_vars, collapse = ", "),
                                       mean_rmse_top_lags = current_rmse))

    remaining <- drop_groupmates(current_vars[length(current_vars)], remaining[-best_i])
    step <- step + 1
  }

  list(selected_vars = current_vars, history = history)
}

## ----------------------------------------------------------------------
## Robustness check: run select_variables_forced() under several
## different FIXED hyperparameter settings
## ----------------------------------------------------------------------
default_gbm_param_settings <- list(
  modest      = list(n.trees = 1000, interaction.depth = 1, shrinkage = 0.01,  n.minobsinnode = 10),
  generous    = list(n.trees = 3000, interaction.depth = 1, shrinkage = 0.005, n.minobsinnode = 10),
  regularized = list(n.trees = 2000, interaction.depth = 1, shrinkage = 0.005, n.minobsinnode = 25)
)

select_variables_robustness_check <- function(sub_data_by_lag, folds, forced_vars, candidate_vars, response_col,
                                               gbm_param_settings = default_gbm_param_settings,
                                               top_n_lags = 5, exclusive_groups = list()) {

  results_by_setting <- names(gbm_param_settings) %>%
    purrr::set_names() %>%
    purrr::map(function(setting_name) {
      select_variables_forced(sub_data_by_lag, folds, forced_vars, candidate_vars, response_col,
                               gbm_params = gbm_param_settings[[setting_name]], top_n_lags = top_n_lags,
                               exclusive_groups = exclusive_groups)
    })

  ## per-candidate-variable stability: how many settings selected/retained it. 
  stability <- candidate_vars %>%
    purrr::map_dfr(function(v) {
      selected_in <- names(results_by_setting) %>%
        purrr::keep(~ v %in% results_by_setting[[.x]]$selected_vars)
      tibble::tibble(variable = v,
                      n_settings_selected = length(selected_in),
                      n_settings_total    = length(results_by_setting),
                      settings_selected_in = paste(selected_in, collapse = ", "))
    }) %>%
    dplyr::arrange(dplyr::desc(n_settings_selected))

  list(results_by_setting = results_by_setting, stability = stability)
}
