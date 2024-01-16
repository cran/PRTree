#' @title Predict method for PRTree
#'
#' @description Predicted values based on a prtree object.
#'
#' @param object Object of class inheriting from \code{"prtree"}
#' @param newdata A matrix with new values for the covariates.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list with the following arguments
#'
#'  \item{yhat}{The predicted values.}
#'
#'  \item{newdata}{The matrix with the covariates new values.}
#'
#' @export
#'
predict.prtree = function(object, newdata, ...){

	mod <- object
	X_test <- newdata

  # variables
  Iindep = 1

  nrow = nrow(X_test)
  ncol = ncol(X_test)
  X_test = as.matrix(X_test)
  inf = mod$regions[, 'inf']
  sup = mod$regions[, 'sup']
  n_terminal_nodes = as.integer(sum(mod$nodes_matrix_info[, 'isTerminal']))
  tn = which(mod$nodes_matrix_info[, 'isTerminal']==1)

  P = matrix(0.0, ncol = n_terminal_nodes, nrow = nrow)
  gamma = mod$gamma
  sigma = as.matrix(mod$sigma)
  yhat_test = numeric(nrow)

  out <- .Fortran("predict_pr_treer", NAOK = TRUE,
                  Iindep = Iindep, nrow = nrow, ncol = ncol,
                  X_test = X_test, inf = inf, sup = sup,
                  n_terminal_nodes = n_terminal_nodes, tn = tn,
                  P = P, gamma = gamma, dim_sigma = dim(sigma),
                  sigma = sigma, yhat_test= yhat_test)


  return(list(yhat = out$yhat,newdata = X_test))
}

