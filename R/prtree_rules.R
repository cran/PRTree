# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. NULL coalescing operator #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.null.default <- function(x, default) {
  if (is.null(x)) default else x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. Rule builder #----
# For parameters with default values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.build <- function(param_name, type, positive = NULL,
                        non_negative = NULL, range = NULL,
                        options = NULL) {
  function(...) {
    list(
      type = type,
      positive = positive,
      non_negative = non_negative,
      range = range,
      length = 1,
      options = options,
      custom = function(x) {
        .null.default(x, .pr.meta[param_name, "default"][[1]])
      }
    )
  }
}
.rule.build.deprecated <- function(param_name, new_param) {
  function(verbose = TRUE, ...) {
    list(
      custom = function(x) {
        if (!is.null(x) && verbose) {
          .warning.message(paste0(
            " Parameter '", param_name, "' is deprecated.\n",
            " Use ", new_param, " instead ."
          ))
        }
        x
      }
    )
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Printing/Saving #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.iprint <- .rule.build(param_name = "iprint", type = "integer")
.rule.verbose <- .rule.build(param_name = "verbose", type = "logical")
.rule.debug <- .rule.build(param_name = "debug", type = "logical")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.y <- function(...) .rule.xy_common(is.X = FALSE)
.rule.X <- function(...) .rule.xy_common(is.X = TRUE)
.rule.xy_common <- function(is.X = TRUE, n_obs = NULL, n_feat = NULL) {
  # set variable name (for printing)
  name <- if (is.X) "`X`" else "`y`"

  # set default messages for type
  allowed_types <- if (is.X) {
    "a matrix, data frame, or numeric vector"
  } else {
    "a numeric vector, a one-column matrix, or a one-column data frame"
  }

  # create a list with a custom validation function for both X and y
  list(
    custom = function(x) {
      # TYPE CHECK - what structures are allowed
      if (!is.matrix(x) && !is.data.frame(x) && !is.vector(x)) {
        .stop.message(paste0(name, " must be ", allowed_types, "."))
      }

      # EMPTY CHECK - no empty objects
      if (length(x) == 0) {
        .stop.message(paste0(name, " cannot be empty."))
      }

      # COLUMN/ROWS COUNT - specific rules for matrix/data.frame
      if (is.matrix(x) || is.data.frame(x)) {
        if (!is.X && ncol(x) > 1) {
          .stop.message(paste0(
            "`y` must have exactly one column when provided as a matrix or data frame."
          ))
        }
      }

      # NA CHECK - y cannot have missing values
      if (!is.X && any(is.na(x))) {
        .stop.message("NA's are not allowed in the response vector 'y'.")
      }

      # NUMERIC CHECK - all data must be numeric
      if (is.data.frame(x)) {
        is_num <- vapply(x, is.numeric, FUN.VALUE = logical(1))
        if (!all(is_num)) {
          non_num <- names(x)[!is_num]
          large <- length(non_num) > 5
          .stop.message(paste0(
            " All columns in ", name, " must be numeric.\n",
            " Categorical variables or factors must be encoded\n",
            " (e.g., one-hot encoding) prior to this step.\n",
            " Non-numeric column(s) detected:",
            if (large) "\n  ",
            paste(non_num, collapse = ", ")
          ))
        }
      } else if (!is.numeric(x)) {
        .stop.message(paste0(name, " must be numeric."))
      }

      # CONVERSION: Convert X to matrix if needed
      #             Ensure y is a flat vector
      if (is.X && !is.matrix(x)) {
        x <- as.matrix(x)
        if (is.null(colnames(x))) {
          colnames(x) <- paste0("X", 1:ncol(x))
        }
      } else if (!is.X) x < as.vector(x)

      # DIMENSION CHECK (for X)
      if (is.X) {
        if (!is.null(n_obs) && nrow(x) != n_obs) {
          .stop.message(sprintf(
            "Dimension mismatch: 'y' has %d observations, but 'X' has %d rows.",
            n_obs, nrow(x)
          ))
        }
        if (!is.null(n_feat) && ncol(x) != n_feat) {
          .stop.message(sprintf(
            "Number of columns in 'newdata' (%d) does not match the training data (%d).",
            ncol(x), n_feat
          ))
        }
      }

      # FORTRAN: set storage mode
      storage.mode(x) <- "double"
      return(x)
    }
  )
}
.rule.n_obs <- function(y = NULL, stratify = NULL, is_NA = NULL, ...) {
  list(
    type = "integer",
    positive = TRUE,
    length = 1,
    custom = function(x) {
      # SAMPLING:
      #   is_NA is required when stratify = TRUE
      #   (validated first and check against n_obs if not NULL)
      if (!is.null(stratify) && stratify) {
        return(length(is_NA))
      }

      # DATA CHECK:
      #   If y exists, n_obs is determined by y
      #   is_NA and y will never be validated togheter.
      if (!is.null(y)) {
        return(length(y))
      }

      # GENERAL CASE:
      #  just check for length 1 + positive + integer.
      x
    }
  )
}
.rule.n_feat <- function(X = NULL, ...) {
  list(
    type = "integer",
    positive = TRUE,
    length = 1,
    custom = function(x) {
      # If X exists, n_feat is determined by X
      if (!is.null(X)) x <- ncol(X)
      x
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grid parameters #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.sigma_grid <- function(n_feat = NULL, ...) {
  list(
    type = "numeric",
    custom = function(x) {
      # NULL is allowed (will be generated later)
      if (is.null(x)) {
        return(x)
      }

      # Validate the object type
      if (!is.vector(x) && !is.matrix(x)) {
        .stop.message(paste0(
          "Invalid format for 'sigma_grid'.\n\n",
          "'sigma_grid' must be one of:\n",
          "   * NULL (to generate grid automatically)\n",
          "   * a numeric vector (for a single sigma value per feature)\n",
          "   * a numeric matrix (for multiple sigma combinations)\n\n",
          "You provided an object of class: ",
          paste(class(x), collapse = ", "), "\n\n",
          "Please correct your input and try again."
        ))
      }

      # Handle matrix case
      if (is.matrix(x)) {
        # If n_feat is available, check compatibility
        if (!is.null(n_feat) && ncol(x) != n_feat) {
          .stop.message(sprintf(
            "The regressor matrix has %d columns, but sigma_grid has %d columns.",
            n_feat, ncol(x)
          ))
        }
        # grid_size_final will be computed later as nrow(x)
        return(x)
      }

      # Handle vector case
      # n_feat might be NULL (not available yet)
      if (is.null(n_feat)) {
        # Just convert to matrix with 1 row
        return(matrix(x, nrow = 1))
      }

      # n_feat is available
      if (length(x) > 1 && length(x) != n_feat) {
        .stop.message(sprintf(
          "The regressor matrix has %d columns, but sigma_grid vector has %d values.\n",
          "Please provide either 1 or %d values.",
          n_feat, length(x), n_feat
        ))
      }

      # Convert to matrix with 1 row and n_feat columns
      matrix(x, nrow = 1, ncol = n_feat)
    }
  )
}
.rule.tiny_sigma <- function(verbose = TRUE, sigma_grid = NULL, ...) {
  list(
    type = "numeric",
    non_negative = TRUE,
    length = 1,
    custom = function(x) {
      # If NULL, just return
      if (is.null(x)) {
        return(x)
      }

      # If sigma_grid is provided, check compatibility
      # Check if first row matches tiny_sigma (within tolerance)
      if (!is.null(sigma_grid)) {
        first_row <- sigma_grid[1, ]
        if (!all(abs(first_row - x) < 1e-15)) {
          if (verbose) {
            .warning.message(sprintf(paste0(
              "'tiny_sigma = %g' does not match the first row of the provided grid.\n",
              "Setting 'tiny_sigma = NULL' for consistency.\n",
              "First row values: %s"),
              x,
              paste(format(first_row, digits = 4), collapse = ", ")
            ))
          }
          return(NULL) # Return NULL when incompatible
        }
        return(x)
      }

      # No sigma_grid - check for large values
      if (x > 1e-10 && verbose) {
        .warning.message(sprintf(paste0(
          "'tiny_sigma = %g' is relatively large.\n",
          "Typical values are very small (e.g., 1e-20)."),
          x
        ))
      }
      x
    }
  )
}
.rule.grid_size <- function(tiny_sigma = NULL, sigma_grid = NULL, verbose = TRUE, ...) {
  list(
    type = "integer",
    non_negative = TRUE,
    length = 1,
    custom = function(x) {
      # If sigma_grid exists, grid_size is irrelevant
      if (!is.null(sigma_grid)) {
        ns <- nrow(sigma_grid)
        return(if(is.null(tiny_sigma)) ns else ns - 1)
      }

      # If NULL, get default from metadata
      x <- .null.default(x, .pr.meta["grid_size", "default"][[1]])

      # grid_size = 0 is only valid with tiny_sigma
      if (x == 0 && is.null(tiny_sigma)) {
        x <- .pr.meta["grid_size", "default"][[1]]
        if (verbose) {
          .warning.message(paste0(
            "'grid_size = 0' is only allowed when 'tiny_sigma' is provided.\n",
            sprintf("Using grid_size = %g (default value).", x)
          ))
        }
      }

      # return number of NORMAL candidates requested
      x
    }
  )
}
.rule.min_mult <- function(sigma_grid = NULL, grid_size, ...) {
  list(
    type = "numeric",
    non_negative = TRUE,
    length = 1,
    custom = function(x) {
      # If sigma_grid exists or grid_size = 1,
      # min_mult is irrelevant - remove from output
      if (!is.null(sigma_grid) || grid_size < 2) {
        return(NULL)
      }

      # Otherwise, apply default if NULL
      .null.default(x, .pr.meta["min_mult", "default"][[1]])
    }
  )
}
.rule.max_mult <- function(min_mult = NULL, ...) {
  list(
    type = "numeric",
    length = 1,
    custom = function(x) {
      # If sigma_grid exists or grid_size = 1,
      # max_mult is set to NULL and min_mult is irrelevant - remove from output
      if (is.null(min_mult)) {
        return(NULL)
      }

      # If NULL, use default
      x <- .null.default(x, .pr.meta["max_mult", "default"][[1]])

      # Validate against min_mult
      if (x <= min_mult) {
        .stop.message("'max_mult' must be greater than 'min_mult'.")
      }
      x
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stopping criteria #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.max_terminal_nodes <- .rule.build(
  param_name = "max_terminal_nodes", type = "integer", positive = TRUE
)
.rule.max_depth <- function(max_terminal_nodes, ...) {
  list(
    type = "integer",
    non_negative = TRUE,
    length = 1,
    custom = function(x) {
      if (is.null(x)) x <- max_terminal_nodes - 1
      x
    }
  )
}
.rule.cp <- .rule.build(param_name = "cp", type = "numeric", non_negative = TRUE)
.rule.n_min <- .rule.build(param_name = "n_min", type = "integer", positive = TRUE)
.rule.perc_x <- .rule.build.deprecated(param_name = "perc_x", new_param = "prop_x")
.rule.prop_x <- function(verbose = TRUE, perc_x = NULL, ...) {
  list(
    type = "numeric",
    range = c(0, 1),
    length = 1,
    custom = function(x) {
      # if none is provided use the default value
      if (is.null(x) && is.null(perc_x)) {
        return(.pr.meta["prop_x", "default"][[1]])
      }

      # case 1: perc_x privided
      #  if prop_x = NULL replace and return
      #  if prop_x is provided warning (keep)
      if (!is.null(perc_x)) {
        if (is.null(x)) {
          if (verbose) {
            cat(sprintf("Note: 'perc_x = %g' provided.\n", perc_x))
            .warning.message(sprintf(paste0(
              " 'perc_test' is deprecated.\n",
              " Setting prop_x =  %g instead."
            ), perc_x))
          }
          return(perc_x)
        } else if (verbose) {
          .warning.message(paste0(sprintf(
            "Both 'perc_x = %g' and 'prop_x = %g' provided.\n", perc_x, x
          ), "'perc_x' is deprecated. Using 'prop_x'."))
        }
      }

      # case 2: only prop_x is provided
      x
    }
  )
}
.rule.p_min <- .rule.build(param_name = "p_min", type = "numeric", range = c(0, 1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Missing data #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.fill_type <- .rule.build(
  param_name = "fill_type", type = "integer",
  options = .pr.meta["fill_type", "options"][[1]]
)
.rule.proxy_crit <- .rule.build(
  param_name = "proxy_crit", type = "character",
  options = .pr.meta["proxy_crit", "options"][[1]]
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Split search #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.n_candidates <- .rule.build(
  param_name = "n_candidates",
  type = "integer", positive = TRUE
)
.rule.by_node <- .rule.build(param_name = "by_node", type = "logical")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribution #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.dist <- .rule.build(
  param_name = "dist", type = "character", options = .pr.meta["dist", "options"][[1]]
)
.rule.dist_pars <- function(dist, ...) {
  all_args <- list(...)
  list(
    type = "list",
    custom = function(x) {
      # Ensure we have a list
      x <- .null.default(x, list())

      # Get parameter info for this distribution
      par_name <- .pr.mapping$dist[[dist]]$par_name

      # If distribution needs a parameter, validate it
      if (!is.null(par_name)) {
        # Priority: direct parameter overrides dist_pars
        if (!is.null(all_args[[par_name]])) {
          x[[par_name]] <- all_args[[par_name]]
        }

        # Get parameter value
        par_value <- x[[par_name]]

        # Validate existence
        if (is.null(par_value)) {
          .stop.message(sprintf(
            "'%s' must be provided for distribution '%s'.\n",
            par_name, dist
          ))
        }

        # Validate value (numeric and positive)
        if (!is.numeric(par_value) || par_value <= 0) {
          .stop.message(sprintf(
            "'%s' must be a positive number.", par_name
          ))
        }
      } else {
        # No parameter needed (e.g., normal)
        par_name <- "par_value"
        par_value <- 0.0
      }

      # Store Fortran info IN THE LIST
      storage.mode(par_value) <- "double"
      x$dist_code <- .pr.mapping$dist[[dist]]$dist_code
      x$par_name <- par_name
      x[[par_name]] <- par_value
      x
    }
  )
}
.rule.param <- function(dist, ...) {
  if (is.null(.pr.mapping$dist[[dist]]$par_name)) {
    return(list())
  }
  list(type = "numeric", positive = TRUE, length = 1)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cross-validation #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.method <- .rule.build(
  param_name = "method", type = "character",
  options = .pr.meta["method", "options"][[1]]
)
.rule.n_rep <- function(method, ...) {
  list(
    type = "integer",
    positive = TRUE,
    length = 1,
    custom = function(x) {
      if (is.null(x)) {
        x <- .pr.meta["n_rep", "default"][[1]]
      } else {
        # validate the number of folds/iterations
        switch(method,
               kfold = if (!is.numeric(x) || length(x) != 1 || x < 2) {
                 .stop.message("k-fold CV: \n'n_rep' (number of folds) must be a single integer >= 2")
               },
               montecarlo = if (!is.numeric(x) || length(x) != 1 || x < 1) {
                 .stop.message("Monte Carlo CV: \n'n_rep' (number of iterations) must be a single integer >= 1")
               }
        )
      }
      x
    }
  )
}
.rule.fold_idx <- function(method, n_rep, n_obs, is.cv = FALSE, ...) {
  list(
    type = "integer",
    custom = function(x) {
      # not provided or not cross validation
      if (is.null(x) || !is.cv) {
        x <- NULL
      } else {
        # if Monte Carlo CV, ignore fold_idx
        if(!method == "kfold") return(NULL)

        # validate the number of folds/iterations
        if(length(unique(x)) != n_rep) {
          .stop.message(" 'fold_idx' must contain `n_rep` unique indexes.")
        }
        if(length(x) != n_obs) {
          .stop.message(" 'fold_idx' must have total length equal to nrow(X).")
        }
      }
      x
    }
  )
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. warning message for cv #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.zero.null.warning <- function(method, type) {
  null <- type == "null"
  .warning.message(paste0(
    " For the current scenario:\n",
    "  * sigma is not fixed, \n",
    "  * `method = ", method, "`,\n",
    "  * `only_sigma = FALSE` \n",
    "  * `prop_valid` = `", ifelse(null, "NULL", "0"), "`\n",
    " The consequences are:\n",
    if (method == "kfold") "  * the reserve fold is used for testing (usual behavior)\n",
    if (null) "  * `prop_valid` will be set to `0`\n",
    "  *  the complete training will be used to select sigma (WARNING)"
  ))
}
.rule.prop_montecarlo <- function(verbose = TRUE, prop_valid = NULL, prop_test = NULL) {
  # Here only_sigma = FALSE:
  # test set size determined by the user
  # validation set determined by the user
  # skip step 1 in fit.tree only when prop_val = 0

  # validate prop_valid
  if (is.null(prop_valid)) {
    zero_valid <- FALSE # (use the default value)
  } else {
    prop_valid <- .check.numeric.range(
      value = prop_valid, range = c(0, 1),
      openlr = c(FALSE, TRUE), name = "prop_valid"
    )
    zero_valid <- !(prop_valid > 0)
    if (zero_valid && verbose) .zero.null.warning(method = "montecarlo", type = "zero")
  }

  # validate prop_test
  if (is.null(prop_test)) {
    zero_test <- FALSE
  } else {
    prop_test <- .check.numeric.range(
      value = prop_test, range = c(0, 1),
      openlr = c(FALSE, TRUE), name = "prop_test"
    )
    zero_test <- !(prop_test > 0)
  }

  # final check
  if (zero_valid && zero_test) {
    .stop.message(paste0(
      " Invalid scenario:\n",
      "  * sigma is not fixed, \n",
      "  * `method = montecarlo`,\n",
      "  * `sigma_only = FALSE` \n",
      " `prop_valid` and `prop_test` cannot be both `NULL` or `0`"
    ))
  }

  # return if at least one is positive!
  return(invisible(TRUE))
}
.rule.prop_kfold <- function(verbose = TRUE, prop_valid = NULL, prop_test = NULL) {
  # Here only_sigma = FALSE:
  # test set size determined by folds
  # validation set determined by the user
  # skip step 1 in fit.tree only when prop_val = 0
  if (is.null(prop_valid)) {
    # if prop_valid was not provided, use default value
    return(invisible(TRUE))
  }

  # validate prop_valid provided by the user
  prop_valid <- .check.numeric.range(
    value = prop_valid, range = c(0, 1),
    openlr = c(FALSE, TRUE), name = "prop_valid"
  )
  if (!prop_valid > 0 && verbose) .zero.null.warning(method = "kfold", type = "zero")
  return(invisible(TRUE))
}
.rule.only_sigma <- function(verbose = TRUE, prop_valid = NULL, prop_test = NULL,
                             sigma_grid = NULL, grid_size, method, ...) {
  list(
    custom = function(x) {
      # CASE 1: sigma is fixed
      #  (there are no candidates to compare or select from)
      only_tiny <- grid_size < 1
      fixed_sigma <- !is.null(sigma_grid) && grid_size < 2

      if (fixed_sigma || only_tiny) {
        # (a) only_sigma = TRUE: not allowed -> convert to only_sigma = FALSE
        # (b) only_sigma = FALSE:
        #      set prop_valid = NULL
        #      validate prop_test based on the method.
        if (verbose && !is.null(x) && x) {
          .warning.message(paste0(
            "  `only_sigma = TRUE` ignored when sigma is fixed.\n",
            "   Setting `only_sigma` to FALSE."
          ))
        }
        # prop_test
        #  k-fold = 1/n_rep (test is done in the fold)
        #  montecarlo = required and validated as a proportion (0,1)
        if (method == "montecarlo" && !is.null(prop_test)) {
          #  if provided, validate as a proportion in (0,1)
          prop_test <- .check.numeric.range(
            value = prop_test, range = c(0, 1),
            openlr = c(TRUE, TRUE), name = "prop_test"
          )

        }
        return(FALSE)
      }

      # If NULL, set to default value
      # If provided, validate as a single logical value
      if (is.null(x)) {
        x <- .pr.meta["only_sigma", "default"][[1]]
      } else if (!is.logical(x) || length(x) != 1) {
        .stop.message("'only_sigma' must be a single logical value")
      }

      # CASE 2: sigma is not fixed and
      #  (a) only_sigma = TRUE:
      #      prop_test = NULL (no test set for any method)
      #      prop_valid
      #       - k-fold = 1/n_rep (test is done in the fold)
      #       - montecarlo = required and validated as a proportion (0,1)
      if (x) {
        if (!is.null(prop_test) && verbose) {
          .warning.message(paste0(
            " `only_sigma = TRUE` overrides `prop_test`.\n",
            " Setting prop_test to NULL (no test set)."
          ))
        }
        if (method == "montecarlo" && !is.null(prop_valid)) {
          if (is.null(prop_valid)) {
            if (verbose) {
              .warning.message(paste0(
                " `only_sigma = FALSE`\n",
                " `method = montecarlo',\n",
                " `prop_valid` cannot be NULL or 0.",
                sprintf(
                  " Using the default value, prop_valid = %g'.",
                  .pr.meta["prop_valid", "default"][[1]]
                )
              ))
            }
          } else { #  if provided, validate as a proportion in (0,1)
            prop_valid <- .check.numeric.range(
              value = prop_valid, range = c(0, 1),
              openlr = c(TRUE, TRUE), name = "prop_valid"
            )
          }
        }
        return(x)
      }

      #  (b) only_sigma = FALSE:
      #      validate prop_valid based on the method.
      #      validate prop_test based on the method.
      rule <- switch(method,
                     kfold = .rule.prop_kfold,
                     montecarlo = .rule.prop_montecarlo
      )
      check <- rule(verbose = verbose, prop_valid = prop_valid, prop_test = prop_test)
      x
    }
  )
}
.rule.stratify <- .rule.build(param_name = "stratify", type = "logical")
.rule.update_final <- .rule.build(param_name = "update_final", type = "logical")

.rule.is_NA <- function(stratify = FALSE, n_obs = NULL, ...) {
  list(
    type = "logical",
    custom = function(x) {
      # CHECK: user provided is_NA?
      has_NA_vec <- !is.null(x)

      # VALIDADE:
      #  stratify = TRUE: compute n_obs from is_NA (required)
      #
      #  stratify = FALSE: check if length(is_NA) = n_obs
      #   - is_NA is used only to compute stistics about missing values.
      #   - n_obs is required
      if (stratify) {
        if (!has_NA_vec) {
          .stop.message("'is_NA' must be provided when stratify = TRUE")
        }
      } else {
        if (is.null(n_obs) || n_obs < 1) .stop.message("'n_obs' must be > 0")
        if (has_NA_vec && (length(x) != n_obs)) {
          .stop.message("When 'is_NA' is provided, its length must equal n_obs")
        }
      }
      x
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data splitting #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.rule.perc_test <- .rule.build.deprecated(param_name = "perc_test", new_param = "prop_hold")
.rule.prop_hold <- function(verbose = TRUE, perc_test = NULL, ...) {
  list(
    length = 1,
    custom = function(x) {
      # if none is provided use the default value
      if (is.null(x) && is.null(perc_test)) {
        return(.pr.meta["prop_hold", "default"][[1]])
      }
      if (!is.null(x)) {
        x <- .check.numeric.range(
          value = x, range = c(0, 1), openlr = c(FALSE, TRUE), name = "prop_hold"
        )
      }

      # case 1: perc_x privided
      #  if prop_x = NULL replace and return
      #  if prop_x is provided warning (keep)
      if (!is.null(perc_test)) {
        perc_test <- .check.numeric.range(
          value = perc_test, range = c(0, 1), openlr = c(FALSE, TRUE),
          name = "perc_test", deprecated = TRUE
        )
        if (is.null(x)) {
          if (verbose) {
            cat(sprintf("Note: 'perc_test = %g' provided.\n", perc_test))
            .warning.message(sprintf(paste0(
              " 'perc_test' is deprecated.\n",
              " Setting prop_hold =  %g instead."
            ), perc_test))
          }
          return(perc_test)
        } else if (verbose) {
          .warning.message(paste0(sprintf(
            "Both 'perc_test = %g' and 'prop_hold = %g' provided.\n", perc_test, x
          ), "'perc_test' is deprecated. Using 'prop_hold'."))
        }
      }

      # case 2: only prop_x is provided
      x
    }
  )
}
.rule.prop_valid <- function(verbose = TRUE, method, sigma_grid = NULL, grid_size,
                             only_sigma, n_rep, ...) {
  # validation is perfomed during the check of only_sigma
  # here values are just updated
  list(
    custom = function(x) {
      # CASE 1: sigma is fixed
      # (there are no candidates to compare or select from)
      only_tiny <- grid_size < 1
      fixed_sigma <- !is.null(sigma_grid) && grid_size < 2
      if (fixed_sigma || only_tiny) {
        return(NULL)
      }

      # CASE 2: only_sigma = TRUE:
      # validate prop_valid based on the method (already done)
      if (only_sigma && method == "kfold") {
        return(1.0 / n_rep)
      }

      # CASE 3: only_sigma = FALSE:
      # validate prop_valid based on the method (already done)
      .null.default(x, .pr.meta["prop_valid", "default"][[1]])
    }
  )
}
.rule.prop_test <- function(verbose = TRUE, method, only_sigma, n_rep, ...) {
  # validation is perfomed during the check of only_sigma
  # here values are just updated
  list(
    custom = function(x) {
      # CASE 1: only_sigma = TRUE
      # (no test set required)
      if (only_sigma) {
        return(NULL)
      }

      # CASE 2: only_sigma = FALSE:
      # validate prop_test based on the method (already done)
      if (method == "kfold") {
        return(1.0 / n_rep)
      }
      .null.default(x, .pr.meta["prop_test", "default"][[1]])
    }
  )
}
.rule.idx_train <- function(n_obs = NULL, is.cv = FALSE, verbose = TRUE, ...) {
  list(
    type = "integer",
    custom = function(x) {
      if (is.null(x)) {
        return(x)
      } else if (is.cv) {
        if (verbose) {
          .warning.message(paste0(
            " `idx_train` provided by ignored.\n",
            " Cross-validation generates its own training/test splits.\n"
          ))
        }
        return(NULL)
      }
      if (!is.null(n_obs) && any(x < 1 | x > n_obs)) {
        .stop.message(" 'idx_train' contains indices outside 1:nrow(X).")
      }
      x
    }
  )
}
