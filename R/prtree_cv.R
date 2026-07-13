# Last revision: July/2026

#' @title Cross-Validation for Probabilistic Regression Trees
#'
#' @description
#' Performs cross-validation for Probabilistic Regression Trees. Supports
#' both k-fold and Monte Carlo (repeated random subsampling) methods.
#'
#' @param y A numeric vector of response values.
#'
#' @param X A numeric matrix or data frame of predictor variables.
#'
#' @param control_cv A list of cross-validation control parameters, typically
#'   created by [pr_tree_control_cv]. Default values are taken from [pr_tree_control_cv]
#'   for any parameters not specified. Alternatively, control parameters can
#'   be passed directly via the `...` argument.
#'
#' @template template_verbose
#' @param ... Control parameters to be passed to [pr_tree_control].
#'
#' @details
#' The [pr_tree_cv] function provides flexible cross-validation for Probabilistic
#' Regression Trees, supporting three common scenarios:
#'
#' **1. Sigma + test error estimation (default)**
#'
#' When both `prop_valid > 0` and a grid of \eqn{\boldsymbol\sigma} values is provided
#' (either via `control$sigma_grid` or automatically generated), the function
#' performs a two-step procedure in each iteration/fold:
#' \itemize{
#'   \item The training data is further split into training and validation sets
#'         according to `prop_valid`.
#'   \item The validation set is used to select the optimal \eqn{\boldsymbol\sigma} from the grid.
#'   \item The model is then refitted using the **full training set** (training +
#'         validation) with the selected \eqn{\boldsymbol\sigma}.
#'   \item Finally, the refitted model is evaluated on the test set.
#' }
#' This yields both the selected \eqn{\boldsymbol\sigma} values for each iteration/fold and
#' estimates of the test error.
#'
#' **2. Test error estimation only (sigma vector fixed)**
#'
#' When `prop_valid = 0` or `NULL`, no internal validation is performed.
#' This is appropriate when \eqn{\boldsymbol\sigma} is already known (e.g., from a previous
#' analysis) and you only need to estimate the test error. In this case:
#' \itemize{
#'   \item The model is trained directly on the full training set using the
#'         provided \eqn{\boldsymbol\sigma} values (a vector or `1` by `n_feat` matrix).
#'   \item The test error is computed on the test set.
#' }
#' This mode is computationally lighter and focuses solely on error estimation.
#'
#' **3. Sigma selection only (no test error)**
#'
#' When `prop_test = 0` (for Monte Carlo) or when you're interested only in the
#' selected \eqn{\boldsymbol\sigma} values and not in test error estimation, the function
#' can be used without a test set. In this case:
#' \itemize{
#'   \item All data is used for training and validation.
#'   \item The optimal \eqn{\boldsymbol\sigma} is selected via internal validation
#'         (if `prop_valid > 0`) or using the full training set
#'         (if `prop_valid = 0`).
#'   \item No test error is computed (`rmse_by_rep` will contain `NULL`).
#' }
#' This mode is useful when you want to estimate \eqn{\boldsymbol\sigma} from the entire
#' dataset before refitting a final model, or when you're only interested in
#' the stability of \eqn{\boldsymbol\sigma} estimates across different splits.
#'
#' **Cross-validation methods**
#'
#' Two resampling methods are available:
#' \itemize{
#'   \item **Monte Carlo CV** (`method = "montecarlo"`): In each of `n_iter`
#'         iterations, the data is randomly split into training (`1 - prop_test`)
#'         and testing (`prop_test`) sets. Observations may appear in multiple
#'         test sets across iterations.
#'   \item **k-fold CV** (`method = "kfold"`): The data is partitioned into
#'         `n_folds` disjoint subsets. Each fold serves as the test set once,
#'         while the remaining `n_folds - 1` folds are used for training and
#'         validation. This ensures every observation is tested exactly once.
#' }
#'
#' @return
#' A list of class `prtree.cv` with components:
#' \item{sigma_matrix}{Matrix of selected \eqn{\boldsymbol\sigma} values,
#'  with one row per iteration/fold and one column per feature.}
#' \item{rmse_by_rep}{A list with components:
#'   \itemize{
#'     \item \code{validation}: Numeric vector of validation RMSE values
#'           (NULL if no validation set used).
#'     \item \code{test}: Numeric vector of test RMSE values
#'           (NULL if `only_sigma = TRUE` or no test set).
#'   }}
#' \item{track_sigma}{Logical indicating whether sigma varies across
#'   iterations/folds (`TRUE`) or is fixed (`FALSE`). When `TRUE`, different
#'   sigma values may be produced for each iteration/fold, either because they
#'   are estimated from the data of that iteration (`sigma_grid = NULL`) or
#'   selected from a grid of candidates (`grid_size > 1`). When `FALSE`, the
#'   same sigma value is used for all iterations/folds.}
#'
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(123)
#' X <- matrix(runif(1000, 0, 10), ncol = 2)
#' y <- 2 * sin(X[, 1]) + 0.5 * X[, 2] + rnorm(500, 0, 0.2)
#'
#' # Default k-fold CV
#' cv1 <- pr_tree_cv(y, X)
#' plot(cv1)
#'
#' # Monte Carlo CV with custom parameters
#' control_cv <- pr_tree_control_cv(
#'   method = "montecarlo",
#'   n_iter = 20,
#'   prop_test = 0.3,
#'   prop_valid = 0.2
#' )
#' cv2 <- pr_tree_cv(y, X, control_cv = control_cv)
#' plot(cv2)
#' }
#'
#' @export
#'
pr_tree_cv <- function(y, X, control_cv = list(), verbose = TRUE,...) {
  # validate the arguments
  .validate.args(y = y, X = X)

  # Ensure y is a flat vector and set to double
  y <- as.vector(y)
  storage.mode(y) <- "double"
  n_obs <- length(y)

  # Convert to matrix if needed
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"

  # Check dimensions
  if (nrow(X) != n_obs) {
    .stop.message(sprintf(
      "Dimension mismatch: 'y' has %d observations, but 'X' has %d rows.",
      n_obs, nrow(X)
    ))
  }

  # Replaces default values by user defined parameters (and validate)
  control_cv <- .update.control(y = y, X = X, control_cv, params = list(...), is.cv = TRUE)

  # Missing values for stratification
  is_NA <- if (control_cv$stratify) apply(is.na(X), 1, any) else NULL

  # Print information
  if (verbose) {
    switch(control_cv$method,
      kfold = .message.kfold(control_cv),
      montecarlo = .message.montecarlo(control_cv)
    )
  }

  # Select internal function based on method
  fun <- switch(control_cv$method, kfold = .cv.kfold, montecarlo = .cv.montecarlo)

  # Run CV
  result <- fun(
    y = y, X = X, control_cv = control_cv, is_NA = is_NA, verbose = verbose
  )

  # add metadata
  result$config <- control_cv[intersect(attr(.pr.meta, "control_cv"), names(control_cv))]

  class(result) <- "prtree.cv"
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Monte Carlo cross-validation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.cv.montecarlo <- function(y, X, is_NA, control_cv, verbose) {
  # get parameters
  n_obs <- control_cv$n_obs
  n_iter <- control_cv$n_rep
  track_sigma <- control_cv$track_sigma
  track_test <- control_cv$track_test

  # Initialize results
  if (track_sigma) {
    # save sigma from each iteration
    sigma_matrix <- matrix(NA, nrow = n_iter, ncol = control_cv$n_feat)
    colnames(sigma_matrix) <- colnames(X)
  } else {
    # sigma is a vector of fixed value
    sigma_matrix <- control_cv$sigma_grid
  }
  sigma_rmse_by_rep <- if (track_sigma) numeric(n_iter) else NULL
  test_rmse_by_rep <- if (track_test) numeric(n_iter) else NULL

  # --- Main loop ---
  for (iter in 1:n_iter) {
    if (verbose) {
      cat(sprintf("\n---------- Iteration %d/%d ----------\n", iter, n_iter))
    }

    # Split the dataset into training and testing sets
    idx_train <- train_test_split(
      n_obs = n_obs,
      prop_test = if (track_test) control_cv$prop_test else 0,
      stratify = control_cv$stratify,
      is_NA = is_NA
    )$idx_train

    if (control_cv$null_sigma) {
      control_cv <- .update.sigma_grid(
        control = control_cv,
        sigma_grid = .expand_sigma_grid(
          X = X,
          n_feat = control_cv$n_feat,
          idx_train = idx_train,
          grid_size = control_cv$grid_size,
          min_mult = control_cv$min_mult,
          max_mult = control_cv$max_mult,
          tiny_sigma = control_cv$tiny_sigma,
          verbose = control_cv$verbose
        )
      )
    }

    # build the tree and select sigma
    # First: build the tree (train) and select sigma (valid)
    # Second: re-build the tree (train + valid) and estimate rmse (test)
    fit <- .fit.tree_rep(
      y = y, X = X, idx_train = idx_train, is_NA = is_NA, control_cv = control_cv
    )

    # Save results
    if (track_sigma) {
      sigma_matrix[iter, ] <- fit$sigma
      sigma_rmse_by_rep[iter] <- sqrt(fit$mse["validation"])
    }
    if (track_test) test_rmse_by_rep[iter] <- sqrt(fit$mse["test"])

    if (verbose) {
      if (track_sigma) {
        cat(sprintf(
          "  Selected sigma: %s\n",
          .format_vector(sigma_matrix[iter, ], max_print = 5, digits = 5)
        ))
        cat(sprintf("  Validation RMSE: %.6f\n", sigma_rmse_by_rep[iter]))
      }
      if (track_test) cat(sprintf("  Test RMSE: %.6f\n", test_rmse_by_rep[iter]))
    }
  }

  list(
    sigma_matrix = sigma_matrix,
    rmse_by_rep = list(
      validation = sigma_rmse_by_rep,
      test = test_rmse_by_rep
    ),
    track_sigma = track_sigma
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: k-fold cross-validation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.cv.kfold <- function(y, X, is_NA, control_cv, verbose) {
  # get parameters
  n_feat <- control_cv$n_feat
  n_folds <- control_cv$n_rep
  track_sigma <- control_cv$track_sigma
  track_test <- control_cv$track_test

  # In case folds were not provided, create folds (stratified if requested)
  if (!is.null(control_cv$fold_idx)) {
    fold_idx <- control_cv$fold_idx
  } else {
    if (control_cv$stratify) {
      fold_idx <- .sample.stratified_folds(control_cv$n_obs, n_folds, is_NA)
    } else {
      fold_idx <- sample(rep(1:n_folds, length.out = control_cv$n_obs))
    }
  }

  # Initialize results
  if (track_sigma) {
    # save sigma from each iteration
    sigma_matrix <- matrix(NA, nrow = n_folds, ncol = n_feat)
    colnames(sigma_matrix) <- colnames(X)
  } else {
    # sigma is a vector of fixed value
    sigma_matrix <- control_cv$sigma_grid
  }
  sigma_rmse_by_fold <- if (track_sigma) numeric(n_folds) else NULL
  test_rmse_by_fold <- if (track_test) numeric(n_folds) else NULL

  for (fold in 1:n_folds) {
    if (verbose) {
      cat(sprintf("\n---------- Fold %d/%d ----------\n", fold, n_folds))
    }

    # Split the dataset into training (k-1 folds) and testing (1 fold) sets
    idx_train <- which(fold_idx != fold)

    if (control_cv$null_sigma) {
      control_cv <- .update.sigma_grid(
        control = control_cv,
        sigma_grid = .expand_sigma_grid(
          X = X,
          n_feat = control_cv$n_feat,
          idx_train = idx_train,
          grid_size = control_cv$grid_size,
          min_mult = control_cv$min_mult,
          max_mult = control_cv$max_mult,
          tiny_sigma = control_cv$tiny_sigma,
          verbose = control_cv$verbose
        )
      )
    }

    # build the tree and select sigma
    # First: build the tree (train) and select sigma (valid)
    # Second: re-build the tree (train + valid) and estimate rmse (test)
    fit <- .fit.tree_rep(
      y = y, X = X, idx_train = idx_train,
      is_NA = is_NA, control_cv = control_cv
    )

    # Save results
    if (track_sigma) {
      sigma_matrix[fold, ] <- fit$sigma
      sigma_rmse_by_fold[fold] <- sqrt(fit$mse["validation"])
    }
    if (track_test) test_rmse_by_fold[fold] <- sqrt(fit$mse["test"])

    if (verbose) {
      if (track_sigma) {
        cat(sprintf(
          "  Selected sigma: %s\n",
          .format_vector(sigma_matrix[fold, ], max_print = 5, digits = 5)
        ))
        cat(sprintf("  Validation RMSE: %.6f\n", sigma_rmse_by_fold[fold]))
      }
      if (track_test) cat(sprintf("  Test RMSE: %.6f\n", test_rmse_by_fold[fold]))
    }
  }

  list(
    sigma_matrix = sigma_matrix,
    rmse_by_rep = list(
      validation = sigma_rmse_by_fold,
      test = test_rmse_by_fold
    ),
    track_sigma = track_sigma
  )
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Initial message for Monte Carlo cross-validation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.message.montecarlo <- function(control_cv) {
  cat("\n========================================================\n")
  cat("PRTree Monte Carlo Cross-Validation\n")
  cat("========================================================\n")
  cat(sprintf("Iterations: %d\n", control_cv$n_rep))

  # Test set information
  if (control_cv$track_test) {
    cat(sprintf("Test set: %.1f%% of data per iteration\n", 100 * control_cv$prop_test))
  } else {
    cat("Test set: none\n")
  }

  # Validation information
  if (!control_cv$track_sigma) {
    # Sigma fixed
    cat("Sigma: fixed\n")
    cat("Validation set: none\n")
  } else if (!control_cv$simplified) {
    # Need sigma + internal validation
    cat(sprintf(
      "Validation set: %.1f%% of training data (internal split)\n",
      100 * control_cv$prop_valid
    ))
    cat("Sigma selection: internal validation split\n")
  } else {
    # Need sigma + without internal validation
    cat("Sigma selection: full training set (no internal validation)\n")
    cat("Validation set: none\n")
  }

  cat(sprintf("Features: %d\n", control_cv$n_feat))
  cat("========================================================\n")
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Initial message for k-fold cross-validation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.message.kfold <- function(control_cv) {
  cat("\n========================================================\n")
  cat("PRTree k-Fold Cross-Validation\n")
  cat("========================================================\n")
  cat(sprintf("Folds: %d\n", control_cv$n_rep))

  # Test set information
  if (control_cv$track_test) {
    cat(sprintf(
      "Test set: 1 fold per iteration (%.1f%% of data)\n",
      100 / control_cv$n_rep
    ))
  } else {
    cat("Test set: none\n")
  }

  # Validation information
  if (!control_cv$track_sigma) {
    # Sigma fixed
    cat("Sigma: fixed\n")
    cat("Validation set: none\n")
  } else if (!control_cv$simplified) {
    # Need sigma + internal validation
    cat(sprintf(
      "Validation set: %.1f%% of training data (internal split)\n",
      100 * control_cv$prop_valid
    ))
    cat("Sigma selection: internal validation split\n")
  } else if (control_cv$only_sigma) {
    # Need sigma + validation in the fold
    cat("Sigma selection: using fold as validation\n")
    cat("Validation set: none\n")
  } else {
    # Need sigma + no internal validation
    cat("Sigma selection: training fold (no internal validation)\n")
    cat("Validation set: none\n")
  }

  cat(sprintf("Features: %d\n", control_cv$n_feat))
  cat("========================================================\n")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Select sigma and re-estimate the tree #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.fit.tree_rep <- function(y, X, idx_train, is_NA, control_cv) {
  n_train <- length(idx_train)
  n_obs <- control_cv$n_obs

  # First call computes the validation MSE and
  # Second call computes the test MSE

  # (a) First call: Train model with combined train + validation data
  #    - Skip if sigma is fixed (no candidates to validate)
  if (control_cv$multiple_sigma) {
    # Split the training set into a validation training set and
    # validation testing set. Update the idx_train
    control_cv <- .update.idx_train(
      control = control_cv,
      new_idx = train_test_split(
        n_obs = n_train,
        prop_test = control_cv$prop_valid,
        stratify = control_cv$stratify,
        is_NA = is_NA[idx_train]
      )$idx_train,
      n_obs = n_train
    )

    # Fit the tree (train set) and compute the mse (validation set)
    # for each sigma in the grid. Returns the one with smaller mse.
    tree <- .pr_tree(
      y = y[idx_train], X = X[idx_train, , drop = FALSE], control = control_cv
    )

    # Save the validation MSE
    if (control_cv$simplified) {
      mse_valid <- unname(tree$mse["train"])
    } else {
      mse_valid <- unname(tree$mse["validation"])
    }

    # Update the sigma_grid to fit the final tree
    control_cv <- .update.sigma_grid(
      control = control_cv,
      sigma_grid = matrix(tree$sigma, nrow = 1)
    )
  } else {
    # update in the second call, if needed
    mse_valid <- NULL
  }

  # update control_pr
  control_cv <- .update.idx_train(control = control_cv, new_idx = idx_train, n_obs = n_obs)

  # (b) Second call:
  #    - use the fitted model and just compute predictions
  #    - fit the model again using all idx_train and predict
  if(control_cv$multiple_sigma && control_cv$track_test && !control_cv$update_final){
    # If the code gets here
    #   (i) the First call was exectued (there are multiple sigmas) to choose from
    #  (ii) there exist a test set
    # (iii) we do not want to update the final tree using all training data.
    #       This, only call predict using the tree fitted in the first call
    idx_test <- setdiff(1:control_cv$n_obs, control_cv$idx_train)
    yhat <- predict.prtree(tree, newdata = X[idx_test, , drop = FALSE])
    tree$mse <- c(tree$mse, test = mean((y[idx_test] - yhat)^2))
  } else {
    # If the code gets here only one of these is TRUE
    #  (i) first call was excecuted and the final model needs update.
    # (ii) there was only one sigma to test
    # In any case, use all training data to fil the model.
    tree <- .pr_tree(y = y, X = X, control = control_cv)
    # Check if the mse_valid needs to be computed.
    # If first call was excecuted, mse_valid was already computed.
    # If first call was not executed (multiple sigmas = FALSE)
    if(!control_cv$multiple_sigma){
      if(control_cv$track_test){
        # no validation set, only train + test
        mse_valid <- unname(tree$mse["train"])
      } else {
        # no test set, only train + validation (named "test" because grid_size = 1)
        mse_valid <- unname(tree$mse["test"])
      }
    }
  }
  return(list(sigma = tree$sigma, mse = c(validation = mse_valid, tree$mse)))
}
