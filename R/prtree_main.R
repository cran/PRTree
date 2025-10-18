#' Probabilistic Regression Trees (PRTrees)
#'
#' @description Fits a Probabilistic Regression Tree (PRTree) model. This is the
#'   main user-facing function of the package.
#'
#' @param y A numeric vector for the dependent variable.
#'
#' @param X A numeric matrix or data frame for the independent variables.
#'
#' @param control A list of control parameters, typically created by
#'   `pr_tree_control()`. Alternatively, control parameters can be passed
#'   directly via the `...` argument.
#'
#' @param ... Control parameters to be passed to `pr_tree_control()`.
#'   These will override any parameters specified in the `control` list.
#'
#' @return An object of class `prtree` containing the fitted model. This is a
#'   list with the following components
#'
#' \item{yhat}{The estimated values for `y`.}
#'
#' \item{XRegion}{A matrix with two columns indicating the terminal node (region)
#'   each observation belongs to. The first column (`TRUE`) may have `NA` for
#'   observations with missing values. The second column (`Internal`) shows the
#'   region assigned by the algorithm.}
#'
#' \item{dist}{The Fortran code corresponding to the distribution used. (For
#' prediction purposes)}
#' \item{par_dist}{Parameters related to the distribution (if any).}
#'
#' \item{fill_type}{Fortran code corresponding to the method used to fill the
#' matrix P when missing values are present.}
#'
#' \item{P}{The matrix of probabilities for each terminal node.}
#'
#' \item{gamma}{The values of the \eqn{\gamma_j} weights estimated for the
#' returned tree}
#'
#' \item{MSE}{The mean squared error for the training, test/validation, and
#'   global datasets.}
#'
#' \item{sigma}{The optimal \eqn{\sigma} vector selected by the grid search.}
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
#' # Fit model with custom controls passed directly
#' reg <- pr_tree(y, X, max_terminal_nodes = 9, perc_test = 0)
#'
#' plot(
#'   X[order(X)], reg$yhat[order(X)],
#'   xlab = "x", ylab = "cos(x)", col = "blue", type = "l"
#' )
#' points(
#'   X[order(X)], y[order(X)],
#'   xlab = "x", ylab = "cos(x)", col = "red"
#' )
#'
#' @importFrom stats sd
#' @importFrom utils modifyList
#'
#' @export
#'
pr_tree <- function(y, X, control = list(), ...) {
  # update the control list
  ctrl <- .update.control(control = control, ...)

  # check for NA in y
  if (any(is.na(y))) stop("NA's not allowed in vector y")
  ctrl$n_obs <- length(y)

  # check the distribution
  ctrl$my_dist <- .get.dist(dist = ctrl$dist, dist_pars = ctrl$dist_pars)

  # check the proxy criterion
  ctrl$crit <- switch(ctrl$proxy_crit,
    mean = 1L,
    var = 2L,
    both = 3L
  )

  # set variables to pass to Fortran
  X <- as.matrix(X)
  ctrl$n_feat <- ncol(X)

  # get the indexes for the training sample
  train <- .get.idx(X = X, idx_train = ctrl$idx_train, perc_test = ctrl$perc_test)
  ctrl[names(train)] <- train

  # Check the sigma_grid argument
  # - if NULL, uses a grid based on the variances of X
  ctrl$sigma_grid <- .get.sigma.grid(X = X, idx_train = ctrl$idx_train, sigma_grid = ctrl$sigma_grid, grid_size = ctrl$grid_size)
  ctrl$grid_size <- nrow(ctrl$sigma_grid)

  # call FORTRAN
  .prtree(X = X, y = y, ctrl = ctrl)
}


.prtree <- function(X, y, ctrl) {
  nc <- 2 * ctrl$max_terminal_nodes - 1
  out <- .Fortran("pr_tree_fort",
    n_obs = ctrl$n_obs,
    n_feat = ctrl$n_feat,
    n_train = ctrl$n_train,
    idx_train = ctrl$idx_train,
    y = y,
    X = X,
    n_sigmas = ctrl$grid_size,
    sigmas = ctrl$sigma_grid,
    int_param = c(
      fill_type = ctrl$fill_type,
      crit = ctrl$crit,
      max_terminal_nodes = ctrl$max_terminal_nodes,
      max_depth = max(1L, ctrl$max_depth),
      n_min = max(1L, ctrl$n_min),
      n_cand = max(1L, ctrl$n_candidates),
      by_node = ifelse(ctrl$by_node, 1L, 0L),
      dist = ctrl$my_dist$dist_code,
      iprint = ctrl$iprint
    ),
    dble_param = c(
      perc_x = ctrl$perc_x,
      p_min = ctrl$p_min,
      cp = ctrl$cp,
      par_dist = ctrl$my_dist$dist_pars
    ),
    n_tn = 1L,
    P = matrix(0, nrow = ctrl$n_obs, ncol = ctrl$max_terminal_nodes),
    gamma = numeric(ctrl$max_terminal_nodes),
    yhat = numeric(ctrl$n_obs),
    MSE = c(training = 0.0, test = 0.0, global = 0.0),
    nodes_matrix_info = matrix(0L, nrow = nc, ncol = 5),
    thresholds = numeric(nc),
    bounds = matrix(0, ncol = 2, nrow = ctrl$n_feat * nc),
    sigma_best = numeric(ctrl$n_feat),
    XRegion = integer(ctrl$n_obs),
    NAOK = TRUE
  )

  # process and return the output
  return(.get.output(object = out, my_dist = ctrl$my_dist))
}

.get.output <- function(object, my_dist) {
  # sample type
  sample <- rep("train", object$n_obs)
  if (object$n_train < object$n_obs) {
    sample[-object$idx_train] <- ifelse(object$n_sigmas == 1, "test", "validation")
  }

  # number of terminal nodes
  n_tn <- object$n_tn
  naux <- 2 * n_tn - 1

  # Probability matrix
  object$P <- object$P[, 1:n_tn, drop = FALSE]
  rownames(object$P) <- sample
  tn <- which(object$nodes_matrix_info[, 2] == 1)
  colnames(object$P) <- paste0("R", 1:n_tn, "(Id:", tn, ")")

  # Format other outputs
  object$gamma <- object$gamma[1:n_tn]
  names(object$yhat) <- sample
  names(object$XRegion) <- sample

  object$nodes_matrix_info <- as.data.frame(
    cbind(
      object$nodes_matrix_info[1:naux, , drop = FALSE],
      threshold = object$thresholds[1:naux]
    )
  )
  colnames(object$nodes_matrix_info) <- c(
    "node", "isTerminal", "fatherNode", "depth", "feature", "threshold"
  )

  regions <- as.data.frame(object$bounds[1:(object$n_feat * naux), , drop = FALSE])
  colnames(regions) <- c("inf", "sup")
  regions$node <- rep(object$nodes_matrix_info$node, each = object$n_feat)
  regions$feature <- rep(1:object$n_feat, naux)
  regions$isTerminal <- rep(object$nodes_matrix_info$isTerminal, each = object$n_feat)
  regions[regions[, "inf"] <= -.Machine$double.xmax, "inf"] <- -Inf
  regions[regions[, "sup"] >= .Machine$double.xmax, "sup"] <- Inf

  object$XRegion <- cbind(object$XRegion, object$XRegion)
  colnames(object$XRegion) <- c("TRUE", "Internal")
  miss <- apply(object$X, 1, function(x) any(is.na(x)))
  if (length(miss) > 0) object$XRegion[miss, "TRUE"] <- NA

  final <- list(
    yhat = object$yhat,
    XRegion = object$XRegion,
    dist = my_dist,
    fill_type = object$int_param["fill_type"],
    P = object$P,
    gamma = object$gamma,
    MSE = object$MSE,
    sigma = object$sigma_best,
    nodes_matrix_info = object$nodes_matrix_info,
    regions = regions[c("node", "feature", "inf", "sup", "isTerminal")]
  )
  class(final) <- c(class(final), "prtree")
  return(final)
}
