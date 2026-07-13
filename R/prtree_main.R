# Last revision: July/2026

#' @title
#' Probabilistic Regression Trees (PRTrees)
#'
#' @description
#' Fits a Probabilistic Regression Tree (PRTree) model. This is the main
#' user-facing function of the package.
#'
#' @details
#' The tree is built using the training data specified by `idx_train` (or
#' determined by `prop_test`). If validation data is available (the complement
#' of `idx_train`), it is used to select the optimal \eqn{\boldsymbol\sigma} values from
#' the grid.
#'
#' **Grid generation**:
#' \itemize{
#'  \item If `sigma_grid` is not provided, it is automatically generated using
#'   **only the training set** (the observations identified by `idx_train`,
#'   which can be passed via `control` or calculated internally from `prop_test`).
#'   The validation/test data is used only for evaluation.
#'  \item Alternatively, users can create a custom grid using the helper function
#'   `expand_sigma_grid(X, ...)` and pass it via the `control` parameter. This
#'   allows using the full dataset (or any desired subset) to define the grid of
#'   \eqn{\boldsymbol\sigma} values to be tested.
#' }
#'
#' **Important**: The returned model is trained **only on the training set**,
#' not on the combined training + validation data. This matches standard
#' practice where validation data is used only for hyperparameter selection, not
#' for final parameter estimation.
#'
#' If you want to retrain the model using the selected \eqn{\boldsymbol\sigma} values on a
#' larger dataset (e.g., combining training and validation), you can:
#' \enumerate{
#'  \item Extract the optimal \eqn{\boldsymbol\sigma} from the fitted model (`model$sigma`)
#'  \item Create a new control object with `sigma_grid = matrix(model$sigma, nrow = 1)`
#'   and set `prop_test = 0` to use all data for training
#'  \item Run [pr_tree] again
#' }
#'
#' @param y A numeric vector for the dependent variable.
#'
#' @param X A strictly numeric matrix or data frame for the independent
#'   variables. Categorical variables or factors must be encoded (e.g., one-hot
#'   encoding) prior to model fitting.
#'
#' @param control A list of control parameters, typically created by
#'   [pr_tree_control]. Default values are taken from [pr_tree_control]
#'   for any parameters not specified. Alternatively, control parameters can
#'   be passed directly via the `...` argument.
#'
#' @param ... Control parameters to be passed to [pr_tree_control].
#'   These will override any parameters specified in the `control` list.
#'
#' @return An object of class `prtree` containing the fitted model. This is a
#'   list with the following components
#'
#' \item{n_obs}{Number of observations used in the model.}
#'
#' \item{n_feat}{Number of features (predictor variables) in the model.}
#'
#' \item{features}{The features names: either `colnames(X)` or a vector with
#' generic names `X1,...,Xn_feat`}
#'
#' \item{yhat}{The estimated values for `y`.}
#'
#' \item{XRegion}{A matrix with two columns indicating the terminal node (region)
#'   each observation belongs to. The first column (`TRUE`) may have `NA` for
#'   observations with missing values. The second column (`Internal`) shows the
#'   region assigned by the algorithm.}
#'
#' \item{dist}{A list with distribution information:
#'   \itemize{
#'     \item \code{dist_name}: Name of the distribution ("norm", "lnorm", "t", "gamma").
#'     \item \code{dist_code}: Integer code used by Fortran (1-4).
#'     \item \code{par_name}: Name of the distribution parameter (e.g., "df", "shape").
#'     \item \code{<par_name>}: Value of the distribution parameter.
#'   }}
#'
#' \item{fill_type}{Fortran code corresponding to the method used to fill the
#' matrix P when missing values are present.}
#'
#' \item{P}{The matrix of probabilities for each terminal node.}
#'
#' \item{gamma}{The values of the \eqn{\gamma_j} weights estimated for the
#' returned tree}
#'
#' \item{mse}{The mean squared error for the training, test/validation, and
#'   global datasets.}
#'
#' \item{sigma}{The optimal \eqn{\boldsymbol\sigma} vector selected by the grid search.}
#'
#' \item{nodes_matrix_info}{A matrix with information for each node of the tree.}
#'
#' \item{regions}{A data frame with the bounds of each variable in each node of
#' the returned tree.}
#'
#' @examples
#' set.seed(1234)
#' X <- matrix(runif(200, 0, 10), ncol = 1)
#' eps <- matrix(rnorm(200, 0, 0.05), ncol = 1)
#' y <- cos(X) + eps
#'
#' # Fit model
#' reg <- pr_tree(y, X, max_terminal_nodes = 9, prop_hold = 0)
#'
#' # Visualize fit (X is 1D, so we can plot the relationship)
#' ord <- order(X)
#' plot(X[ord], reg$yhat[ord],
#'   type = "l", col = "blue", lwd = 2,
#'   xlab = "x", ylab = "y",
#'   main = "PRTree Fit to cos(x)"
#' )
#' points(X[ord], y[ord], col = "red", cex = 0.5)
#' legend("topright",
#'   legend = c("Fitted", "Observed"),
#'   col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 1)
#' )
#'
#' # Diagnostic plots
#' par(mar = c(3, 4, 1.5, 1))
#' plot(reg, y, which = c(1, 4, 2, 5), ncol = 2)
#' 
#' # Plotting the final tree
#' plot_tree(reg, max_bins = 10)
#'
#' @seealso
#' [pr_tree_cv] for cross-validation options.
#'
#' [pr_tree_control] for setting control parameters.
#'
#' [expand_sigma_grid] for creating custom sigma grids.
#'
#' @export
#'
pr_tree <- function(y, X, control = list(), ...) {
  # validate the data (this will also convert X to the correct format if needed)
  .validate.args(y = y, X = X)

  # update the control list (convert to required format)
  control <- .update.control(y = y, X = X, control = control, params = list(...))

  # call FORTRAN
  .pr_tree(X = X, y = y, control = control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Call the Fortran subroutine #----
# and process and returns the processed output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.pr_tree <- function(X, y, control) {
  nc <- 2 * control$max_terminal_nodes - 1
  keep_full <- control$output_mode == 1
  if(control$n_obs > control$n_train){
    idx_test <- setdiff(1:control$n_obs, control$idx_train)
  } else {
    idx_test = 1
  }
  n_test = control$n_obs - control$n_train

  out <- .Fortran(
    "pr_tree_fort",
		# --------- INPUTS ---------
    n_obs = control$n_obs,
    n_feat = control$n_feat,
    n_train = control$n_train,
    y_train = y[control$idx_train],
	  y_test = y[idx_test],
    X_train = X[control$idx_train, , drop = FALSE],
	  X_test = X[idx_test, , drop = FALSE],
    n_sigmas = control$grid_size_final,
    sigmas = control$sigma_grid,
    int_param = c(
      fill_type = control$fill_type,
      crit = .pr.mapping$proxy_crit[[control$proxy_crit]],
      max_terminal_nodes = control$max_terminal_nodes,
      max_depth = control$max_depth,
      n_min = control$n_min,
      n_cand = control$n_candidates,
      by_node = ifelse(control$by_node, 1L, 0L),
      dist = control$dist_pars$dist_code,
      iprint = control$iprint,
      output_mode = control$output_mode
    ),
    dble_param = c(
      prop_x = control$prop_x,
      p_min = control$p_min,
      cp = control$cp,
      par_dist = control$dist_pars[[control$dist_pars$par_name]]
    ),
    n_tn = 1L,
    # --------- OUTPUTS ---------
    P_train = matrix(0.0,
               nrow = if (keep_full) control$n_train else 1,
               ncol = if (keep_full) control$max_terminal_nodes else 1),
    P_test = matrix(0.0,
               nrow = if (keep_full) max(1, n_test) else 1,
               ncol = if (keep_full) control$max_terminal_nodes else 1),
    gamma = numeric(if (keep_full) control$max_terminal_nodes else 1),
    yhat_train = numeric(if (keep_full) control$n_train else 1),
    yhat_test = numeric(if (keep_full) max(1, n_test) else 1),
    mse = c(train = 0.0, test = 0.0, global = 0.0),
    nodes_matrix_info = matrix(0L,
                              nrow = if (keep_full) nc else 1,
                              ncol = if (keep_full) 5 else 1),
    thresholds = numeric(if (keep_full) nc else 1),
    sigma_best = numeric(control$n_feat),
    XRegion_train = integer(if (keep_full) control$n_train else 1),
    NAOK = TRUE
  )

  # process and return the output
  control$dist_pars$dist_name <- control$dist
  out$idx_test <- idx_test
  out$idx_train <- control$idx_train
  return(.get.output.pr_tree(object = out, dist = control$dist_pars, debug = control$debug))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Processes the output from Fortran #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.output.pr_tree <- function(object, dist, debug) {
  sigma_fixed <- object$n_sigmas == 1
  names(object$mse)[2] <- ifelse(sigma_fixed, "test", "validation")

  if (object$int_param["output_mode"] == 0) {
  return(list(
    mse = object$mse,
    sigma = object$sigma_best
  ))
 }

  # sample type
  sample <- rep("train", object$n_obs)
  if (object$n_train < object$n_obs) {
    sample[-object$idx_train] <- ifelse(sigma_fixed, "test", "validation")
  }

  # number of terminal nodes
  n_tn <- object$n_tn
  naux <- 2 * n_tn - 1

  # Probability matrix
  if(object$n_obs > object$n_train){
    P <- matrix(0, nrow = object$n_obs, ncol = n_tn)
    P[object$idx_train,] <- object$P_train[, 1:n_tn, drop = FALSE]
    P[object$idx_test,] <- object$P_test[, 1:n_tn, drop = FALSE]
    object$P <- P
    object$yhat <- numeric(object$n_obs)
    object$yhat[object$idx_train] <- object$yhat_train
    object$yhat[object$idx_test] <- object$yhat_test
    rm(P)
  } else {
    object$P <- object$P_train[, 1:n_tn, drop = FALSE]
    object$yhat <- object$yhat_train
  }
  object$P_train <- NULL
  object$P_test <- NULL
  object$yhat_train <- NULL
  object$yhat_test <- NULL

  rownames(object$P) <- sample
  tn <- which(object$nodes_matrix_info[, 2] == 1)
  colnames(object$P) <- paste0("R", 1:n_tn, "(Id:", tn, ")")

  # Format other outputs
  object$gamma <- object$gamma[1:n_tn]
  names(object$yhat) <- sample

  object$nodes_matrix_info <- as.data.frame(
    cbind(
      object$nodes_matrix_info[1:naux, , drop = FALSE],
      threshold = object$thresholds[1:naux]
    )
  )
  colnames(object$nodes_matrix_info) <- c(
    "node", "isTerminal", "fatherNode", "depth", "feature", "threshold"
  )

  # Reconstruct regions without needing the bounds array from Fortran
  inf_bounds <- matrix(-Inf, nrow = naux, ncol = object$n_feat)
  sup_bounds <- matrix(Inf, nrow = naux, ncol = object$n_feat)

  if (naux > 1) {
    for (i in 2:naux) {
      father <- object$nodes_matrix_info$fatherNode[i]
      split_feat <- object$nodes_matrix_info$feature[father]
      thr <- object$nodes_matrix_info$threshold[father]

      inf_bounds[i, ] <- inf_bounds[father, ]
      sup_bounds[i, ] <- sup_bounds[father, ]

      # Even indices are left children, odd indices are right children
      if (i %% 2 == 0) {
        sup_bounds[i, split_feat] <- thr
      } else {
        inf_bounds[i, split_feat] <- thr
      }
    }
  }

  regions <- data.frame(
    node = rep(object$nodes_matrix_info$node, each = object$n_feat),
    feature = rep(1:object$n_feat, naux),
    inf = as.vector(t(inf_bounds)),
    sup = as.vector(t(sup_bounds)),
    isTerminal = rep(object$nodes_matrix_info$isTerminal, each = object$n_feat)
  )

  object$XRegion <- matrix(0, nrow = object$n_obs, ncol = 2)
  object$XRegion[object$idx_train, 1] <- object$XRegion_train
  object$XRegion[object$idx_train, 2] <- object$XRegion_train
  colnames(object$XRegion) <- c("TRUE", "Internal")
  rownames(object$XRegion) <- sample
  miss <- logical(object$n_obs)
  miss[object$idx_train] <- apply(object$X_train, 1, function(x) any(is.na(x)))
  if(object$n_obs > object$n_train){
    miss[object$idx_test] <- apply(object$X_test, 1, function(x) any(is.na(x)))
  }
  if (any(miss)) object$XRegion[miss, "TRUE"] <- NA

  features <- colnames(object$X_train)
  if(is.null(features)) features <- paste0("X", 1:ncol(object$X_train))
    
  final <- list(
    n_obs = object$n_obs,
    n_feat = object$n_feat,
    features = colnames(object$X_train),
    yhat = object$yhat,
    XRegion = object$XRegion,
    dist = dist,
    fill_type = object$int_param["fill_type"],
    P = object$P,
    gamma = object$gamma,
    mse = object$mse,
    sigma = object$sigma_best,
    nodes_matrix_info = object$nodes_matrix_info,
    regions = regions[c("node", "feature", "inf", "sup", "isTerminal")]
  )
  if (.null.default(debug, FALSE)) {
    attr(final, "control") <- object[[attr(.pr.meta, "control_pr")]]
  }
  class(final) <- "prtree"
  return(final)
}
