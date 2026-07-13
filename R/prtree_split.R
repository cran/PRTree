# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generates train/test indexes #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Split data into training and testing sets
#'
#' @description Creates indices for training and testing sets. Can use simple
#'   random sampling or stratified sampling based on missing values.
#'
#' @param n_obs Integer. Total number of observations. Required if `stratify =
#'   FALSE`. Ignored when `stratify = TRUE` (in which case `n_obs` is determined
#'   from `is_NA`).
#'
#' @param prop_test Numeric between 0 and 1. Proportion of data to use for
#'   testing. Default is 0.2.
#'
#' @param stratify Logical. If `TRUE`, performs stratified sampling based on
#'   missing values. In this case, `is_NA` must be provided. Default is `FALSE`.
#'
#' @param is_NA Optional logical vector. `TRUE` indicates rows with missing
#'   values. Required if `stratify = TRUE`. If provided when `stratify = FALSE`,
#'   its length must equal `n_obs` and missing statistics will be computed (see
#'   details).
#'
#' @param n_rep Optional integer. Number of folds for creating a
#'   stratified or simple random \eqn{k}-fold partition. If `NULL` (default),
#'   the function generates a single training/test split according to
#'   `prop_test`. When specified, `prop_test` is ignored and the function
#'   returns a vector of fold assignments (`fold_idx`) instead of training
#'   and testing indices.
#'   
#' @details The function creates indices for training and testing sets using
#'   either simple random sampling or stratified sampling based on missing
#'   values.
#'
#' **Stratified sampling:**
#'   When `stratify = TRUE`, the sampling preserves the proportion of rows with
#'   missing values in both training and testing sets. In this case:
#' \itemize{
#'   \item `is_NA` must be provided
#'   \item `n_obs` is ignored (determined from `is_NA`)
#' }
#'
#' **Missing values statistics:**
#'   If `is_NA` is provided but `stratify = FALSE`, the function will still
#'   compute missing values statistics for the resulting splits. The length of
#'   `is_NA` is against `n_obs` to ensure consistency. These statistics will be
#'   available in the [summary] output, with a note indicating that they were
#'   not used for stratification.
#'
#' **Reproducibility:** To ensure reproducible splits, call [set.seed] before
#'   using this function.
#'
#' @return An object of class `idx.split`.
#'
#' If `n_rep = NULL` (default), the returned object contains:
#' \itemize{
#'   \item `idx_train`: Integer vector of training indices.
#'   \item `idx_test`: Integer vector of testing indices (or `NULL` if
#'   `prop_test = 0`).
#'   \item `n_train`: Integer. Number of training observations.
#'   \item `n_test`: Integer. Number of testing observations.
#'   \item `stratified`: Logical indicating whether stratified sampling was
#'   used.
#'   \item `missing_stats`: List with missing-value counts (if `is_NA` is
#'   provided), containing:
#'   \itemize{
#'     \item `train_na`: Number of rows with missing values in the training set.
#'     \item `test_na`: Number of rows with missing values in the testing set.
#'     \item `total_na`: Total number of rows with missing values.
#'   }
#' }
#'
#' If `n_rep` is specified, the returned object contains:
#' \itemize{
#'   \item `fold_idx`: Integer vector indicating the fold assignment of each
#'   observation.
#'   \item `n_rep`: Integer. Number of folds.
#'   \item `stratified`: Logical indicating whether stratified sampling was
#'   used.
#'   \item `fold_stats`: Data frame summarizing each fold, including the number
#'   of observations, the number of rows with missing values, and the
#'   corresponding proportion of missing values.
#' }
#'
#' @examples
#' \donttest{
#' # Simple random split
#' set.seed(123)
#' split <- train_test_split(n_obs = 100, prop_test = 0.3)
#' str(split)
#'
#' # Stratified split by missing values
#' missing <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.2, 0.8))
#' split <- train_test_split(
#'   prop_test = 0.3,
#'   stratify = TRUE,
#'   is_NA = missing
#' )
#' summary(split)
#'
#' # Compute missing stats without stratification
#' missing2 <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.2, 0.8))
#' split <- train_test_split(
#'   n_obs = 100,
#'   prop_test = 0.3,
#'   stratify = FALSE,
#'   is_NA = missing2
#' )
#' summary(split) # Note about missing stats will appear
#' }
#' 
#' # Stratified 10-fold partition
#' folds <- train_test_split(
#'   n_obs = 1000,
#'   is_NA = sample(c(TRUE, FALSE), prob = c(0.1, 0.9), replace = TRUE, size = 1000),
#'   stratify = TRUE,
#'   n_rep = 5
#' )
#' summary(folds)
#' @export
train_test_split <- function(n_obs = NULL, prop_test = 0.2, stratify = FALSE,
                             is_NA = NULL, n_rep = NULL) {
  # validate parameters (output to the local environment)
  prop_hold <- prop_test
  
  if (!is.null(n_rep)) .validate.args(n_rep = n_rep, method = "kfold")
  .validate.args(
    n_obs = n_obs, prop_hold = prop_hold,
    stratify = stratify, is_NA = is_NA
  )
  
  # number of folds not provided.
  # split the indexes into training and test
  if (is.null(n_rep)) {
    # split the indexes
    result <- .sample.idx(
      idx_train = NULL, is_NA = is_NA, n_obs = n_obs,
      prop_test = prop_hold, stratify = stratify
    )
    result$idx_test <- setdiff(1:n_obs, result$idx_train)
    result$n_test <- n_obs - result$n_train
    result$stratified <- stratify
    
    # Add missing stats if provided
    if (!is.null(is_NA)) {
      result$missing_stats <- list(
        train_na = as.integer(sum(is_NA[result$idx_train])),
        test_na = if (result$n_test > 0) as.integer(sum(is_NA[result$idx_test])) else 0L,
        total_na = as.integer(sum(is_NA))
      )
    }
    class(result) <- "idx.split"
    return(result)
  } 
  
  # number of folds provided,
  # create the folds
  if (stratify) {
    fold_idx <- .sample.stratified_folds(n_obs = n_obs, n_folds = n_rep, is_NA = is_NA)
  } else {
    fold_idx <- sample(rep(1:n_rep, length.out = n_obs))
  }
  
  result <- list(
    fold_idx = as.integer(fold_idx),
    n_rep = as.integer(n_rep),
    stratified = stratify
  )
  
  # Add missing statistics if provided
  if (!is.null(is_NA)) {
    n_na <- tapply(is_NA, fold_idx, sum)
    n_obs_fold <- tabulate(fold_idx, nbins = n_rep)
    result$fold_stats <- data.frame(
      fold = seq_len(n_rep),
      n_obs = as.integer(n_obs_fold),
      n_na = as.integer(n_na),
      prop_na = n_na / n_obs_fold
    )
  }
  
  class(result) <- "idx.split"
  return(result)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Process the indexes for the training set#----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.sort.idx <- function(idx_train) {
  idx_train <- as.integer(sort(idx_train))
  list(
    idx_train = idx_train,
    n_train = length(idx_train)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: performs stratified samping #----
# to preserve missing value proportions, as per documentation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.sample.stratified <- function(is_NA, prop_test) {
  # Sampling from the complete group
  idx_complete <- which(!is_NA)
  n_complete <- length(idx_complete)
  n_train_complete <- round(n_complete * (1 - prop_test))

  if (n_train_complete > 1) {
    train_idx_complete <- sample(idx_complete, n_train_complete)
  } else if (n_train_complete == 1) {
    train_idx_complete <- idx_complete
  } else {
    train_idx_complete <- NULL
    n_train_complete <- 0
  }

  # Sampling from the group with NA
  n_NA <- length(is_NA) - n_complete
  n_train_NA <- round(n_NA * (1 - prop_test))

  if (n_train_NA > 1) {
    idx_NA <- which(is_NA)
    train_idx_NA <- sample(idx_NA, n_train_NA)
  } else if (n_train_NA == 1) {
    idx_NA <- which(is_NA)
  } else {
    train_idx_NA <- NULL
    n_train_NA <- 0
  }

  list(
    n_train = as.integer(n_train_complete + n_train_NA),
    idx_train = as.integer(sort(c(train_idx_complete, train_idx_NA)))
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Extract/process the indexes for the training set #----
# if the index for the training sample is provided, ignore prop_test
# just convert to the required format and return
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.sample.idx <- function(idx_train = NULL, is_NA = NULL, n_obs = 0,
                        prop_test = 0, stratify = TRUE) {
  # check if indexes were provided
  if (!is.null(idx_train)) {
    return(.sort.idx(idx_train))
  }

  # If prop_test is provided, use prop_test to obtain the training sample
  if (prop_test > 0 && prop_test < 1) {
    # if stratify, perform stratified random sampling
    if (stratify) {
      return(.sample.stratified(is_NA, prop_test))
    }
    # simple random sampling
    n_train <- max(floor(n_obs * (1 - prop_test)), 1)
    return(.sort.idx(sample(1:n_obs, size = n_train)))
  }

  # Use all data for training if prop_test is 0
  list(
    n_train = as.integer(n_obs),
    idx_train = as.integer(1:n_obs)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Create stratified k-fold indices #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.sample.stratified_folds <- function(n_obs, n_folds, is_NA) {
  # Edge case: no stratification needed
  # If is_NA is NULL, all values are NA, or none are NA,
  # fall back to simple random sampling
  if (is.null(is_NA) || all(is_NA) || !any(is_NA)) {
    return(sample(rep(1:n_folds, length.out = n_obs)))
  }

  # Separate indices by missing status
  idx_complete <- which(!is_NA)
  n_complete <- length(idx_complete)
  idx_missing <- which(is_NA)
  n_missing <- length(idx_missing)

  # Create fold assignments for each group separately
  # rep(1:n_folds, length.out = n) creates a balanced distribution
  # where each fold gets either floor(n/n_folds) or ceiling(n/n_folds)
  # observations, with the remainder distributed to the first folds
  folds_complete <- rep(1:n_folds, length.out = n_complete)
  folds_missing <- rep(1:n_folds, length.out = n_missing)

  # Create the fold index vector and assign values
  # Shuffle within each group to ensure random ordering before fold assignment
  fold_idx <- integer(n_obs)
  fold_idx[idx_complete] <- sample(folds_complete)
  fold_idx[idx_missing] <- sample(folds_missing)

  return(fold_idx)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Print method for indexes #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Print method for idx.split objects
#'
#' @description
#' Prints a brief summary of a train/test split result.
#'
#' @param x An object of class `idx.split`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
print.idx.split <- function(x, ...) {
  
  ## K-fold partition
  if (!is.null(x$fold_idx)) {
    
    cat("K-fold partition:\n")
    cat(sprintf("  Number of folds: %d\n", x$n_rep))
    
    if (!is.null(x$stratified) && x$stratified) {
      cat("  Stratified by missing values\n")
    }
    
    if (!is.null(x$fold_stats)) {
      cat("\n")
      print(x$fold_stats, row.names = FALSE)
    }
    
    return(invisible(x))
  }
  
  ## Holdout partition
  n_obs <- x$n_train + x$n_test
  
  cat("Train/test split:\n")
  cat(sprintf(
    "  Training: %d observations (%.1f%%)\n",
    x$n_train, 100 * x$n_train / n_obs
  ))
  
  if (x$n_test > 0) {
    cat(sprintf(
      "  Testing:  %d observations (%.1f%%)\n",
      x$n_test, 100 * x$n_test / n_obs
    ))
  } else {
    cat("  Testing:  none\n")
  }
  
  if (!is.null(x$stratified) && x$stratified) {
    cat("  Stratified by missing values\n")
  }
  
  invisible(x)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summary method for indexes #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Summary method for `idx.split` objects
#'
#' @description
#' Provides a detailed summary of a train/test split result, including
#' statistics about missing values if available.
#'
#' @param object An object of class `idx.split`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of class `summary.idx.split` containing detailed statistics.
#'
#' @export
summary.idx.split <- function(object, ...) {
  # K-fold partition
  if (!is.null(object$fold_idx)) {
    
    result <- list(
      n_obs = length(object$fold_idx),
      n_rep = object$n_rep,
      stratified = object$stratified,
      fold_stats = object$fold_stats
    )
    
    class(result) <- "summary.idx.split"
    return(result)
  }
  
  # Holdout partition
  n_train <- object$n_train
  n_test <- object$n_test
  n_obs <- object$n_train + n_test

  result <- list(
    n_obs = n_obs,
    n_train = n_train,
    n_test = n_test,
    prop_train = n_train / n_obs,
    prop_test = if (n_test > 0) n_test / n_obs else 0,
    stratified = object$stratified,
    missing_stats = object$missing_stats
  )

  if (!is.null(object$missing_stats)) {
    result$missing_props <- list(
      train = if (n_train > 0) object$missing_stats$train_na / n_train else 0,
      test = if (n_test > 0) object$missing_stats$test_na / n_test else 0,
      total = object$missing_stats$total_na / n_obs
    )
  }

  class(result) <- "summary.idx.split"
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Print summary method for indexes #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Print method for summary of `idx.split` objects
#'
#' @description
#' Prints a detailed summary of a train/test split result.
#'
#' @param x An object of class `summary.idx.split`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
print.summary.idx.split <- function(x, ...) {
  cat("\n========================================================\n")
  
  # ------------------------------------------------------------------
  # K-fold partition
  # ------------------------------------------------------------------
  if (!is.null(x$fold_stats)) {
    
    cat("K-fold Partition Summary\n")
    cat("========================================================\n")
    
    cat(sprintf("Sampling method: %s\n",
      ifelse(x$stratified, "Stratified (by missing values)", "Simple random")))
    
    cat(sprintf("Total observations: %d\n", x$n_obs))
    cat(sprintf("Number of folds: %d\n", x$n_rep))
    
    if (!is.null(x$fold_stats)) {
      cat("\n--- Fold Statistics ---\n")
      
      tab <- x$fold_stats
      tab$prop_na <- sprintf("%.1f%%", 100 * tab$prop_na)
      names(tab) <- c("Fold", "Size", "Missing", "Missing (%)")
      print(tab, row.names = FALSE, right = TRUE)
    }
    
    cat("========================================================\n")
    
    return(invisible(x))
  }
  
  # ------------------------------------------------------------------
  # Holdout partition
  # ------------------------------------------------------------------
  cat("Train/Test Split Summary\n")
  cat("========================================================\n")

  # Sampling method
  cat(sprintf(
    "Sampling method: %s\n",
    ifelse(x$stratified, "Stratified (by missing values)", "Simple random")
  ))

  # Note if missing stats are shown but not used for stratification
  if (!is.null(x$missing_stats) && !x$stratified) {
    cat("\nNote: Missing values statistics are shown below but were NOT\n")
    cat("      used for stratification (stratify = FALSE).\n\n")
  }

  # Basic counts
  cat(sprintf("Total observations: %d\n", x$n_obs))
  cat(sprintf("Training set: %d (%.1f%%)\n", x$n_train, 100 * x$prop_train))

  if (x$n_test > 0) {
    cat(sprintf("Testing set: %d (%.1f%%)\n", x$n_test, 100 * x$prop_test))
  } else {
    cat("Testing set: NULL (no test set)\n")
  }

  # Missing values statistics
  if (!is.null(x$missing_stats)) {
    cat("\n--- Missing Values ---\n")
    cat(sprintf(
      "Training set: %d (%.1f%%)\n",
      x$missing_stats$train_na, 100 * x$missing_props$train
    ))

    if (x$n_test > 0) {
      cat(sprintf(
        "Testing set:  %d (%.1f%%)\n",
        x$missing_stats$test_na, 100 * x$missing_props$test
      ))
    }

    cat(sprintf(
      "Overall:      %d (%.1f%%)\n",
      x$missing_stats$total_na, 100 * x$missing_props$total
    ))
  }

  cat("========================================================\n")
  return(invisible(x))
}
