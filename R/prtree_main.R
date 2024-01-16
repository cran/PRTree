#' Probabilistic Regression Trees (PRTrees)
#'
#' @param y a numeric vector corresponding to the dependent variable
#' @param X A numeric vector, matrix or dataframe corresponding to the independent variables, with the same number of observations as \code{y}.
#' @param sigma_grid optionally, a numeric vector with candidate values for the parameter \eqn{\sigma}, to be passed to the grid search algorithm. If \code{NULL}, the standard deviations of the columns in \code{X} are used. The default is \code{NULL}.
#' @param max_terminal_nodes a non-negative integer. The maximum number of regions in the output tree. The default is 15.
#' @param cp a positive numeric value. The complexity parameter. Any split that does not decrease the MSE by a factor of \code{cp} will be ignored. The default is 0.01.
#' @param max_depth a non-negative integer. The maximum depth of the decision tree. The depth is defined as the length of the longest path from the root to a leaf. The default is 5.
#' @param n_min a non-negative integer, The minimum number of observations in a final node. The default is 5.
#' @param perc_x a positive numeric value. Given a column of \eqn{P}, the value \code{perc_x} is the percentage of rows in this column that must have a probability higher than the threshold \code{p_min} for a splitting attempt to be made in the corresponding region. The default is 0.1.
#' @param p_min a positive numeric value. A threshold probability that controls the splitting process. A splitting attempt is made in a given region only when the proportion of rows with probability higher than \code{p_min}, in the corresponding column of the matrix \eqn{P}, is equal to \code{perc_x}. The default is 0.05.
#'
#' @return
#'
#' \item{yhat}{the estimated values for \code{y}}
#'
#' \item{P}{the matrix of probabilities calculated with the observations in \code{X} for the returned tree}
#'
#' \item{gamma}{the values of the \eqn{\gamma_j} weights estimated for the returned tree}
#'
#' \item{MSE}{the mean squared error calculated for the returned tree}
#'
#' \item{sigma}{the \eqn{\sigma} of the returned tree}
#'
#' \item{nodes_matrix_info}{information related to each node of the returned tree}
#'
#' \item{regions}{information related to each region of the returned tree}
#'
#' @examples
#'
#'set.seed(1234)
#'X = matrix(runif(200, 0, 10), ncol = 1)
#'eps = matrix(rnorm(200, 0, 0.05), ncol = 1)
#'y =  matrix(cos(X) + eps, ncol = 1)
#'reg = PRTree::pr_tree(y, X, max_terminal_nodes = 9)
#'plot(X[order(X)], reg$yhat[order(X)], xlab = 'x', ylab = 'cos(x)', col = 'blue', type = 'l')
#'points(X[order(X)], y[order(X)], xlab = 'x', ylab = 'cos(x)', col = 'red')
#'
#' @export
#'
pr_tree = function(y, X,
                   sigma_grid = NULL,
                   max_terminal_nodes = 15L,
                   cp = 0.01,
                   max_depth = 5L,
                   n_min = 5L,
                   perc_x = 0.1,
                   p_min = 0.05
){

  # variables: setting the maximum size to pass to the Fortran subroutine.

  nrow = nrow(X)
  ncol = ncol(X)
  P = matrix(0, nrow = nrow, ncol = max_terminal_nodes)
  gamma = numeric(max_terminal_nodes)
  yhat = numeric(length(y))
  nodes_matrix_info = matrix(0L ,ncol = 5, nrow = -1 + 2*max_terminal_nodes)
  cutpoints = numeric(-1 + 2*max_terminal_nodes)
  inf = numeric(ncol*(-1 + 2*max_terminal_nodes))
  sup = numeric(ncol*(-1 + 2*max_terminal_nodes))
  XRegion = integer(nrow)

  # making sure X is a matrix
  X = as.matrix(X)

  # if NULL, uses a grid based on the variances of X
  if(is.null(sigma_grid)){
    ss <- apply(X, 2, function(x){stats::sd(x, na.rm = TRUE)})
    sigma <- matrix(seq(min(ss, na.rm = TRUE), max(ss, na.rm = TRUE), length = 5), ncol = 1)
  } else{sigma = matrix(sigma_grid, ncol = 1)}

  out <- .Fortran("pr_treer", NAOK = TRUE, y = y, X = X,
                  nrow = length(y), ncol = ncol(X),
                  sigma = sigma, dim_sigma = dim(sigma),
                  max_terminal_nodes = as.integer(max_terminal_nodes),
                  cp = cp, max_depth = as.integer(max_depth),
                  n_min = as.integer(n_min),
                  perc_x = perc_x, p_min = p_min,
                  Iindep = 1L,  # for now this is just a dummy argument
                  P = P, dim_P = dim(P),
                  gamma = gamma, yhat = yhat, MSE = 0.0,
                  nodes_matrix_info = nodes_matrix_info,
                  cutpoints = cutpoints, inf = inf, sup = sup,
                  sigma_best = 0.0, XRegion = XRegion)

  out$P <- out$P[1:out$dim_P[1], 1:out$dim_P[2]]
  out$gamma <- out$gamma[1:out$dim_P[2]]

  out$nodes_matrix_info = cbind(matrix(out$nodes_matrix_info[1:(-1+2*out$dim_P[2]), ], ncol = 5),
                                out$cutpoints[1:(-1+2*out$dim_P[2])])
  colnames(out$nodes_matrix_info) = c('node', 'isTerminal', 'fatherNode', 'depth', 'varCut', 'cutpoints')

  regions =
    cbind(
      rep(out$nodes_matrix_info[,'node'], each = ncol),
      rep(1:ncol, -1+2*out$dim_P[2]),
      out$inf[1:(ncol*(-1+2*out$dim_P[2]))],
      out$sup[1:(ncol*(-1+2*out$dim_P[2]))],
      rep(out$nodes_matrix_info[,'isTerminal'], each = ncol)
    )
  colnames(regions) = c('node', 'var', 'inf', 'sup', 'isTerminal')
  regions[regions[,'inf'] <= .Machine$double.xmin, 'inf'] = -Inf
  regions[regions[,'sup'] >= .Machine$double.xmax, 'sup'] = Inf

  final = list(yhat = out$yhat, P = out$P, gamma = out$gamma, MSE = out$MSE,
               sigma = out$sigma_best,
               nodes_matrix_info = out$nodes_matrix_info,
               regions = regions)
  class(final) = c(class(final), "prtree")
  return(final)
}
