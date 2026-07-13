# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generates a control list #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Control parameters for PRTree
#'
#' @description This function creates a list of control parameters for the
#'   `pr_tree` function, with validation for each parameter.
#'
#' @template template_sigma_custom
#' @template template_sigma_auto
#' @template template_tree_params
#' @template template_data_split
#' @template template_algo_params
#' @template template_verbose
#'
#' @return A list of class `prtree.control` containing the validated control
#'   parameters.
#'
#' @examples
#' # Get default control parameters
#' controls <- pr_tree_control()
#'
#' # Customize some parameters
#' ctrl1 <- pr_tree_control(max_depth = 5, n_candidates = 5)
#'
#' # equivalent calls
#' ctrl2.v1 <- pr_tree_control(dist = "t", df = 4)
#' ctrl2.v2 <- pr_tree_control(dist = "t", dist_pars = list(df = 4))
#'
#' @export
pr_tree_control <- function(
    sigma_grid = NULL, grid_size = 8, min_mult = 0, max_mult = 2, tiny_sigma = NULL,
    max_terminal_nodes = 15L, max_depth = max_terminal_nodes - 1, cp = 0.01, n_min = 5L,
    prop_x = 0.1, p_min = 0.05, prop_hold = 0.2, idx_train = NULL, fill_type = 2L,
    proxy_crit = "both", n_candidates = 3L, by_node = FALSE, dist = "norm",
    iprint = -1, verbose = TRUE, ...) {
  # --- Parameter Validation ---
  control <- .pr_tree_control(control = .get.input(...))
  class(control) <- "prtree.control"
  return(control)
}

#' @title Print method for PRTree control objects
#'
#' @description
#' Prints a summary of the control parameters for PRTree model fitting.
#'
#' @param x An object of class `prtree.control`, as returned by [pr_tree_control].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
#'
print.prtree.control <- function(x, ...) {
  cat("\n========================================================\n")
  cat("PRTree Control Parameters\n")
  cat("========================================================\n")

  # Sigma grid information
  if (is.null(x$sigma_grid)) {
    cat(sprintf("Sigma grid: not provided (will be generated with size %d)\n", x$grid_size))
    cat(sprintf("Multiplier range: [%g, %g]\n", x$min_mult, x$max_mult))
  } else {
    cat(sprintf("Sigma grid: user-provided (%d candidate vectors)\n", nrow(x$sigma_grid)))
  }
  # Tiny sigma information
  if (!is.null(x$tiny_sigma)) {
    cat(sprintf("Tiny sigma candidate: %g (added as extra row)\n", x$tiny_sigma))
  }

  # Hold-out proportion - role depends on grid_size
  # If null, skip printing (to be used with CV)
  if (!is.null(x$prop_hold)) {
    if (x$prop_hold > 0) {
      grid_size_final <- .grid_size_final(
        sigma_grid = x$sigma_grid,
        grid_size = x$grid_size,
        tiny_sigma = x$tiny_sigma
      )
      if (grid_size_final > 1) {
        role <- "Validation"
      } else {
        role <- "Test"
      }
      cat(sprintf("%s set proportion: %.2f\n", role, x$prop_hold))
    } else {
      cat("Hold-out set: none\n")
    }
  }

  # Training indices info
  if (!is.null(x$idx_train)) {
    cat(sprintf("Training indices: user-provided (%d observations)\n", length(x$idx_train)))
  }

  # Other parameters
  cat(sprintf("Max terminal nodes: %d\n", x$max_terminal_nodes))
  cat(sprintf("Max depth: %d\n", x$max_depth))
  cat(sprintf("Complexity parameter (cp): %.3f\n", x$cp))
  cat(sprintf("Minimum observations per node: %d\n", x$n_min))
  cat(sprintf("Proportion threshold (prop_x): %.2f\n", x$prop_x))
  cat(sprintf("Probability threshold (p_min): %.3f\n", x$p_min))
  cat(sprintf("Fill type: %d\n", x$fill_type))
  cat(sprintf("Proxy criterion: %s\n", x$proxy_crit))
  cat(sprintf("Number of candidates: %d\n", x$n_candidates))
  cat(sprintf("By node: %s\n", toupper(as.character(x$by_node))))
  cat(sprintf("Distribution: %s\n", x$dist))
  if (x$dist_pars$par_name != "par_value") {
    cat(sprintf(
      "  with %s = %.4f\n",
      x$dist_pars$par_name,
      x$dist_pars[[x$dist_pars$par_name]]
    ))
  }

  cat("========================================================\n")
  return(invisible(x))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generates a control list (CV) #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Control parameters for PRTree with optional cross-validation
#'
#' @description Creates a list of control parameters for Probabilistic
#' Regression Trees. The base function [pr_tree_control] returns parameters for
#' tree construction. The extended function [pr_tree_control_cv] returns both
#' tree parameters and cross-validation settings.
#'
#' @details \strong{Choosing the right configuration:}
#'
#' \tabular{lll}{
#' \strong{Goal} \tab \strong{Method} \tab \strong{Parameters} \cr
#' Estimate sigma + test error \tab any \tab `prop_valid > 0`, `only_sigma = FALSE` \cr
#' Test error only (sigma known) \tab any \tab `prop_valid = 0`, `only_sigma = FALSE` \cr
#' Sigma selection only \tab montecarlo \tab `prop_test = 0`, `only_sigma = TRUE` \cr
#' Sigma selection only \tab kfold \tab `only_sigma = TRUE` \cr
#' }
#'
#' @template template_cv_params
#'
#' @param ... Control parameters for PRTree. See [pr_tree_control]
#'
#' @return A list of class `prtree.control_cv` containing merged with
#'   cross-validation settings
#'
#' @seealso
#' [pr_tree] for fitting a single tree,
#'
#' [pr_tree_cv] for cross-validation.
#'
#' @examples
#' # Default k-fold CV
#' cv1 <- pr_tree_control_cv()
#'
#' # Monte Carlo CV with custom parameters
#' cv2 <- pr_tree_control_cv(
#'   method = "montecarlo",
#'   n_iter = 20,
#'   prop_test = 0.3,
#'   prop_valid = 0.2,
#'   stratify = TRUE
#' )
#'
#' @export
#'
pr_tree_control_cv <- function(
    method = "montecarlo", n_rep = 10, only_sigma = FALSE,
    prop_test = 0.2, prop_valid = 0.2, stratify = FALSE,
    update_final = TRUE, fold_idx = NULL,...) {
  # validate parameters
  control <- .pr_tree_control(control = .get.input(...), is.cv = TRUE)
  class(control) <- "prtree.control_cv"
  return(control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Validate the control argument. #----
# Set defaults and convert to Fortran format
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom stats setNames
.pr_tree_control <- function(control, is.cv = FALSE) {
  # All parameters to validate
  dots <- attr(.pr.meta, "dots_remove")

  params <- if (is.cv) "control_cv" else "control_pr"

  all_params <- unique(c(dots, attr(.pr.meta, params)))

  # Get validation order from metadata
  check_order <- attr(.pr.meta, "check_order")
  params_in_order <- intersect(check_order, all_params)

  # Validate ALL parameters in order - rules handle everything!
  ctrl <- do.call(.validate.args, c(
    setNames(control[params_in_order], params_in_order),
    list(is.cv = is.cv, envir = NULL)
  ))

  # Convert environment to list for return
  # remove deprecated ones
  ctrl[dots] <- NULL

  # Set storage modes using metadata attributes
  ctrl[attr(.pr.meta, "pr_int")] <- lapply(
    ctrl[attr(.pr.meta, "pr_int")], .set.storage_mode, "integer"
  )
  ctrl[attr(.pr.meta, "pr_dble")] <- lapply(
    ctrl[attr(.pr.meta, "pr_dble")], .set.storage_mode, "double"
  )

  # return formatted list
  return(ctrl)
}

#' @title Print method for PRTree control objects
#'
#' @description
#' Prints a summary of the control parameters for PRTree model fitting.
#'
#' @param x An object of class `prtree.control_cv`, as returned by
#' [pr_tree_control_cv].
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the original object.
#'
#' @export
#'
print.prtree.control_cv <- function(x, ...) {
  # Print tree parameters (using base control print method)
  # remove prop_test to avoid printing twice
  class(x) <- "prtree.control"
  print(x)
  class(x) <- "prtree.control_cv"

  # Print CV-specific parameters
  cat("Cross-Validation Settings\n")
  cat("========================================================\n")
  cat(sprintf("Method: %s\n", x$method))

  if (x$method == "kfold") {
    cat(sprintf("Number of folds: %d\n", x$n_rep))
    if (x$only_sigma) {
      cat("  Test set: none (sigma selection only).\n")
      cat(sprintf("  Validation (for sigma): 1 fold per iteration (~%.1f%% of data).\n", 100 / x$n_rep))
    } else {
      cat(sprintf("  Test set: 1 fold per iteration (~%.1f%% of data).\n", 100 / x$n_rep))
      if (!is.null(x$prop_valid)) {
        if (x$prop_valid > 0) {
          cat(sprintf("  Validation (for sigma): %.2f of the training folds.\n", x$prop_valid))
        } else {
          cat("  Validation (for sigma): none (sigma selected on training folds).\n")
        }
      } else {
        cat("  Validation (for sigma): none.\n")
      }
    }
  } else { # montecarlo
    cat(sprintf("Number of iterations: %d\n", x$n_rep))
    if (x$only_sigma) {
      cat("  Test set: none (sigma selection only)\n")
    } else if (!is.null(x$prop_test) && x$prop_test > 0) {
      cat(sprintf("  Test proportion: %.2f\n", x$prop_test))
    } else {
      cat("  Test set: none\n")
    }
    if (!is.null(x$prop_valid)) {
      if (x$prop_valid > 0) {
        cat(sprintf("  Validation proportion: %.2f (internal split)\n", x$prop_valid))
      } else {
        cat("  Validation proportion: none (sigma selected on training data)\n")
      }
    } else {
      cat("  Validation: none\n")
    }
  }

  cat(sprintf("Stratification: %s\n", toupper(as.character(x$stratify))))
  cat("========================================================\n")
  return(invisible(x))
}
