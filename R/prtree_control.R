#' Set Control Parameters for PRTree
#'
#' This function creates a list of control parameters for the `pr_tree`
#' function, with validation for each parameter.
#'
#' @param sigma_grid Optional, a numeric value, vector or a matrix with
#'   candidate values for the parameter \eqn{\boldsymbol\sigma}, to be passed to
#'   the grid search algorithm.  If a single numeric value is provided, the code
#'   assumes that \eqn{\sigma_j  = \sigma} for all \eqn{j} covariates. If
#'   \code{NULL}, the standard deviations of the columns in \code{X} are used to
#'   create a grid with values in the interval \eqn{(0, 2\hat\sigma_j]}, with
#'   increments of \eqn{\hat\sigma_j/4}, where \eqn{\hat\sigma_j} denotes the
#'   standard deviation of the \eqn{j}th covariate.  The default is \code{NULL}.
#'
#' @param grid_size Optional, the number of candidate values for `sigma` to
#'   generate when `sigma_grid` is `NULL`. Default is 8.
#'
#' @param max_terminal_nodes A non-negative integer. The maximum number of
#'   regions in the output tree. The default is 15.
#'
#' @param cp A positive numeric value. The complexity parameter. Any split that
#'   does not decrease the MSE by a factor of `cp` will be ignored. The
#'   default is 0.01.
#'
#' @param max_depth A non-negative integer. The maximum depth of the decision
#'   tree. The depth is defined as the length of the longest path from the root
#'   to a leaf. The default is 14.
#'
#' @param n_min A positive integer, The minimum number of observations in a
#'   final node. The default is `max_terminal_nodes - 1`.
#'
#' @param perc_x A positive numeric value between 0 and 1. Given any
#'   column of \eqn{P}, `perc_x` is the minimum proportion of rows that
#'   must have a probability higher than `p_min` for a
#'   splitting attempt to be made in the corresponding region. The split will be
#'   ignored if any of the resulting regions do not meet the same criterion. The
#'   default is 0.1.
#'
#' @param p_min A positive numeric value. A threshold probability that controls
#'   the splitting process. A splitting attempt is made in a given region only
#'   when the proportion of rows with probability higher than `p_min`, in
#'   the corresponding column of the matrix \eqn{P}, is equal to \code{perc_x}.
#'   The default is 0.05.
#'
#' @param perc_test A numeric value between 0 (inclusive) and 1 (exclusive) that
#'   specifies the proportion of the data to be held out for model validation or
#'   testing. Default is 0.2. The role of this hold-out set depends on the
#'   `sigma_grid`
#'   \itemize{
#'     \item **Validation Set:** If `sigma_grid` contains multiple candidate
#'       \eqn{\sigma} values (`grid_size > 1`), `perc_test` of the data is
#'       used as a validation set to select the best \eqn{\sigma} based on
#'       out-of-sample Mean Squared Error (MSE). If `perc_test` is 0, `sigma`
#'       will be selected based on the MSE for the training sample.
#'
#'     \item **Test Set:** If a single, fixed \eqn{\sigma} is provided
#'       (`grid_size = 1`), `perc_test` of the data is used as a test set to
#'       evaluate the final model's performance. If `perc_test` is 0, the
#'       entire dataset is used for training.
#'   }
#'   The data split is performed using stratified sampling to ensure that the
#'   proportion of observations with missing values is similar across the
#'   training and validation/test sets.
#' @param idx_train Indexes for the training sample. Default is `NULL`, in which
#'   case the indexes are computed based on the `perc_test` argument. If
#'   `idx_train` is provided, `perc_test` is ignored.
#'
#' @param fill_type Integer indicating the method to be used to fill the
#'   probability matrix when `X` has NA's. Default is 2.
#'   \itemize{
#'    \item `0`: uniform (same probability for both child nodes).
#'    \item `1`: attributes all probability to the child node that is compatible
#'    with the observed values.
#'    \item `2`: computes the probability restricted to the observed entries
#'   }
#'
#' @param proxy_crit Character. Default is `"both"`. Criterion used to associate
#'   an observation with missing values to a region:
#'   \itemize{
#'    \item `"mean"`: maximizes the difference in means after a split.
#'    \item `"var"`: maximizes the variability between nodes.
#'    \item `"both"`: combines the `"mean"` and `"var"` criteria.
#'   }
#'
#' @param n_candidates Integer. The number of competing candidates to consider
#'   when searching for the best split. To select the candidates, a proxy
#'   improvement measure is used. Then a full analysis is performed to choose
#'   the best among the `n_candidates` candidates. Default is 3.
#'
#' @param by_node Logical. If `TRUE`, the algorithm selects `n_candidates` for
#'   each node and then makes a full analysis to choose the best among all
#'   nodes. Otherwise the `n_candidates` are selected globally. Default is
#'   `FALSE`.
#'
#' @param dist Character. The distribution to be used in the model. One of
#'   `"norm"` (Gaussian), `"lnorm"` (log-normal), `"t"` (Student's \eqn{t}),
#'   or `"gamma"` (Gamma). Default is `"norm"`.
#'
#' @param iprint Integer. Controls the verbosity of the Fortran backend.
#'   Default is -1 (silent).
#'   \itemize{
#'    \item `iprint < 0`: No printing.
#'    \item `iprint = 0`: Prints basic information.
#'    \item `iprint > 0`: As for `iprint = 0` plus progress reports.
#'     }
#'
#' @param ... Extra parameters to be passed to the chosen distribution.
#' \itemize{
#'  \item `"norm"`: Uses the standard Gaussian distribution. No extra
#'   parameters required.
#'  \item `"lnorm"`: Uses the log-normal distribution with `meanlog = 0`.
#'   Requires `sdlog`.
#'  \item `"t"`: Uses the \eqn{t} distribution. Requires `df`.
#'  \item  `"gamma"`: Uses the gamma distribution with `scale = 1`. Requires
#'   `shape`.
#' }
#'
#'
#' @return A list of class `prtree.control` containing the validated control parameters.
#'
#' @export
#' @examples
#' # Get default control parameters
#' controls <- pr_tree_control()
#'
#' # Customize some parameters
#' custom_controls <- pr_tree_control(max_depth = 5, n_candidates = 5)
#'
pr_tree_control <- function(sigma_grid = NULL,
                            grid_size = 8,
                            max_terminal_nodes = 15L,
                            cp = 0.01,
                            max_depth = max_terminal_nodes - 1,
                            n_min = 5L,
                            perc_x = 0.1,
                            p_min = 0.05,
                            perc_test = 0.2,
                            idx_train = NULL,
                            fill_type = 2L,
                            proxy_crit = "both",
                            n_candidates = 3L,
                            by_node = FALSE,
                            dist = "norm",
                            iprint = -1, ...) {
  # --- Parameter Validation ---
  if (!is.null(sigma_grid) && !is.numeric(sigma_grid)) stop("'sigma_grid' must be numeric or NULL.")
  if (grid_size <= 0) stop("'grid_size' must be positive.")
  if (max_terminal_nodes < 1) stop("'max_terminal_nodes' must be at least 1.")
  if (cp < 0) stop("'cp' must be non-negative.")
  if (max_depth < 0) stop("'max_depth' must be non-negative.")
  if (n_min < 1) stop("'n_min' must be at least 1.")
  if (perc_x < 0 || perc_x > 1) stop("'perc_x' must be in [0, 1].")
  if (p_min < 0) stop("'p_min' must be non-negative.")
  if (perc_test < 0 || perc_test >= 1) stop("'perc_test' must be in [0, 1).")
  if (!is.null(idx_train) && !is.numeric(idx_train)) stop("'idx_train' must be a numeric vector or NULL.")
  if (!(fill_type %in% 0:2)) stop("'fill_type' must be 0, 1, or 2.")
  if (!(proxy_crit %in% c("mean", "var", "both"))) stop("'proxy_crit' must be 'mean', 'var', or 'both'.")
  if (n_candidates < 1) stop("'n_candidates' must be at least 1.")
  if (!is.logical(by_node)) stop("'by_node' must be logical.")
  if (!(dist %in% c("norm", "lnorm", "t", "gamma"))) {
    stop(paste0(
      "Distribution '", dist, "' not implemented yet.",
      "\nAvailable choices are:",
      "\n * 'norm' (Gaussian)",
      "\n * 'lnorm' (log-normal)",
      "\n * 't' (Student's t)",
      "\n * 'gamma' (Gamma)"
    ))
  }

  # --- Collect parameters into a list ---
  control_list <- list(
    sigma_grid = sigma_grid,
    grid_size = as.integer(grid_size),
    max_terminal_nodes = as.integer(max_terminal_nodes),
    cp = cp,
    max_depth = as.integer(max_depth),
    n_min = as.integer(n_min),
    perc_x = perc_x,
    p_min = p_min,
    perc_test = perc_test,
    idx_train = idx_train,
    fill_type = as.integer(fill_type),
    proxy_crit = proxy_crit,
    n_candidates = as.integer(n_candidates),
    by_node = by_node,
    dist = dist,
    iprint = as.integer(iprint),
    dist_pars = list(...)
  )

  class(control_list) <- "prtree.control"
  return(control_list)
}



.update.control <- function(control, ...) {
  # Get parameters in ...
  all_params <- list(...)

  # Separate control parameters from distribution parameters
  control_params <- names(formals(pr_tree_control))
  control_args <- all_params[names(all_params) %in% control_params]
  dist_args <- all_params[!names(all_params) %in% c(control_params, "dist_pars")]

  # Merge with explicit control list
  ctrl <- modifyList(control, control_args)

  # Re-validate the ctrl object
  ctrl <- do.call(pr_tree_control, ctrl)

  # copy distribution related parameters (if any)
  pars <- control$dist_pars
  dist_pars <- if (is.null(pars)) list() else pars
  pars <- all_params$dist_pars
  dist_pars <- if (is.null(pars)) dist_pars else modifyList(dist_pars, pars)
  ctrl$dist_pars <- modifyList(dist_pars, dist_args)

  return(ctrl)
}
