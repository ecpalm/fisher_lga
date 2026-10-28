
###############################################################################
# 
# Script title: GBM fitting tools
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
## Tools for fitting GBMs used in 03_fit_straight_line.R
## (iteration 0) and 04_lcp_iteration_cycle.R (iterations 1-5).
## ============================================================================

# Load required packages
sapply(c("gbm", "caret"),
       require, character.only = T)

default_gbm_grid <- expand.grid(
  interaction.depth = c(1, 2),
  n.trees            = seq(100, 3000, 10),
  shrinkage          = 0.005,
  n.minobsinnode     = c(5, 10, 15, 25, 50)
)

## sub_data must be row-aligned with folds$index/indexOut 
fit_gbm_cv <- function(sub_data, predictor_cols, response_col, folds,
                        t_grid = default_gbm_grid, seed = 1234, bag_fraction = 0.5) {

  ## Drops tuneGrid rows that are infeasible for this subgroup's SMALLEST fold. 
  ## Without dropping these rows, performance will be NA. 
  min_fold_train_n <- min(lengths(folds$index))
  feasible <- t_grid$n.minobsinnode <= floor(bag_fraction * min_fold_train_n / 2^t_grid$interaction.depth)
  if (!all(feasible)) {
    message("fit_gbm_cv(): dropping ", sum(!feasible), " of ", nrow(t_grid),
            " tuneGrid row(s) infeasible for this subgroup's smallest fold (",
            min_fold_train_n, " training rows).")
    t_grid <- t_grid[feasible, , drop = FALSE]
  }

  n_folds <- folds$n_folds
  seeds <- vector(mode = "list", length = n_folds + 1)
  set.seed(seed)
  for (i in seq_len(n_folds)) seeds[[i]] <- rep(seed, nrow(t_grid))
  seeds[[n_folds + 1]] <- seed

  ctrl <- caret::trainControl(
    method          = "cv",
    index           = folds$index,
    indexOut        = folds$indexOut,
    allowParallel   = TRUE,
    savePredictions = "final",
    returnResamp    = "none",
    seeds           = seeds
  )

  caret::train(
    x = sub_data[, predictor_cols],
    y = sub_data[[response_col]],
    method    = "gbm",
    trControl = ctrl,
    metric    = "RMSE",
    tuneGrid  = t_grid,
    verbose   = FALSE
  )
}
