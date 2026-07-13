# Last revision: July/2026

#' @title Predict from a Probabilistic Regression Tree Model
#'
#' @description Obtains predictions from a fitted `prtree` object.
#'
#' @param object An object of class `prtree`, as returned by [pr_tree].
#'
#' @param newdata A numeric matrix or data frame containing new data for
#'   which to generate predictions. Must contain the exact same number and
#'   order of predictor variables as the data used to fit the model.
#'   It can contain missing values (`NA`).
#'
#' @param complete Logical. If `FALSE` (default), only the vector of predicted
#'   values is returned. If `TRUE`, a list containing both the predicted values
#'   and the probability matrix `P` is returned.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return If `complete = FALSE`, a numeric vector of predicted values (`yhat`).
#'   If `complete = TRUE`, a list containing:
#'   \item{yhat}{The numeric vector of predicted values.}
#'   \item{P}{The probability matrix for the new data.}
#'
#' @export
#'
predict.prtree <- function(object, newdata, complete = FALSE, ...) {
  # validate categorical/character columns and convert to the correct format
  newdata <- .check.X(newdata, n_feat = object$n_feat)
  nr <- nrow(newdata)
  tn <- which(object$nodes_matrix_info[, "isTerminal"] == 1)
  n_terminal_nodes <- length(tn)

  if (n_terminal_nodes == 1) {
    out <- list(
      yhat = rep(object$gamma, nr),
      P = matrix(1.0, ncol = 1, nrow = nr)
    )
  } else {
    # important features will be extracted directly in fortran
    out <- .Fortran("predict_pr_tree_fort",
      dist = object$dist$dist_code,
      par_dist = object$dist[[object$dist$par_name]],
      fill_type = object$fill_type,
      n_obs = nr,
      n_feat = object$n_feat,
      X_test = newdata,
      thresholds = object$nodes_matrix_info$threshold,
      n_terminal_nodes = n_terminal_nodes,
      tn = tn,
      nodes_info = apply(
        object$nodes_matrix_info[, c("isTerminal", "fatherNode", "depth", "feature")],
        2, as.integer),
      P = matrix(0.0, nrow = nr, ncol = n_terminal_nodes),
      gamma = object$gamma,
      sigma = object$sigma,
      yhat = numeric(nr),
      NAOK = TRUE
    )
  }
  if (!complete) {
    return(out$yhat)
  }
  return(list(yhat = out$yhat, P = out$P))
}
