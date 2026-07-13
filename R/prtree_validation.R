# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: check value against range #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.numeric.range <- function(value, range, openlr, name, deprecated = FALSE) {
  # Check if argument is numeric
  if (!is.numeric(value)) {
    # If an argument is deprecated and has a non-numeric value,
    # return NULL to signal that it should be ignored, preventing further
    # validation errors for a parameter that is no longer in use.
    if (deprecated) {
      return(NULL)
    }
    # if not deprecated, throw an error about type (not range)
    if (!is.null(range)) {
      .stop.message(sprintf("'%s' must be numeric to check range.", name))
    } else {
      # if no range is provided, just check if it's numeric
      return(value)
    }
  }
  # If no range is provided, just return the value (after confirming it's numeric)
  if (is.null(range)) {
    return(value)
  }
  # Check if value is within the specified range
  if (is.null(openlr)) {
    if (any(value < range[1] | value > range[2])) {
      .stop.message(sprintf(
        "'%s' must be in range [%g, %g].", name, range[1], range[2]
      ))
    }
  } else {
    # Open/closed based on openlr
    left_ok <- if (openlr[1]) all(value > range[1]) else all(value >= range[1])
    right_ok <- if (openlr[2]) all(value < range[2]) else all(value <= range[2])
    if (!left_ok || !right_ok) {
      left_br <- if (openlr[1]) "(" else "["
      right_br <- if (openlr[2]) ")" else "]"
      .stop.message(sprintf(
        "'%s' must be in %s%g, %g%s.",
        name, left_br, range[1], range[2], right_br
      ))
    }
  }
  value
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Validate a single argument argument #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.validate.arg <- function(value, name, type = NULL, length = NULL,
                          range = NULL, openlr = NULL,
                          options = NULL, positive = NULL,
                          non_negative = NULL, required = NULL,
                          custom = NULL, msg = NULL) {
  # CASE 1: value is NULL
  if (is.null(value)) {
    # If there's custom validation, let it handle the NULL
    if (!is.null(custom)) {
      if (!is.function(custom)) {
        .stop.message("Custom validation must be a function.")
      }
      value <- custom(value)
    }

    # custom may transform NULL into something
    # If still NULL after custom, fall through to required check
    if (is.null(value) && !is.null(required) && required) {
      .stop.message(sprintf("'%s' is required but was NULL.", name))
    }
    # Not required, return NULL
    return(value)
  }

  # Check type
  if (!is.null(type)) {
    type_ok <- switch(type,
      "numeric" = is.numeric(value),
      "integer" = is.numeric(value) && all(value == as.integer(value)),
      "character" = is.character(value),
      "logical" = is.logical(value),
      "factor" = is.factor(value),
      "data.frame" = is.data.frame(value),
      "matrix" = is.matrix(value),
      "list" = is.list(value),
      "function" = is.function(value),
      FALSE
    )
    if (!type_ok) {
      .stop.message(sprintf("'%s' must be of type '%s'.", name, type))
    }
  }

  # Check length
  if (!is.null(length)) {
    if (length(value) != length) {
      .stop.message(sprintf("'%s' must have length %d.", name, length))
    }
  }

  # Check numeric range
  value <- .check.numeric.range(value = value, range = range, openlr = openlr, name = name)

  # Check positive/non-negative
  if (!is.null(positive) && positive) {
    if (!is.numeric(value)) {
      .stop.message(sprintf("'%s' must be numeric to check positivity.", name))
    }
    if (any(value <= 0)) {
      .stop.message(sprintf("'%s' must be positive.", name))
    }
  } else if (!is.null(non_negative) && non_negative) {
    if (!is.numeric(value)) {
      .stop.message(sprintf("'%s' must be numeric to check non-negativity.", name))
    }
    if (any(value < 0)) {
      .stop.message(sprintf("'%s' must be non-negative.", name))
    }
  }

  # Check allowed options
  if (!is.null(options)) {
    if (!all(value %in% options)) {
      .stop.message(sprintf(
        "'%s' must be one of: %s.", name, paste(options, collapse = ", ")
      ))
    }
  }

  # Custom validation
  if (!is.null(custom)) {
    if (!is.function(custom)) {
      .stop.message("Custom validation must be a function.")
    }
    value <- custom(value)
  }

  return(value)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. Validate multiple arguments at once #----
# (stops at first error)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.validate.args <- function(is.cv = FALSE, envir = parent.frame(),
                           rules.mapping = attr(.pr.meta, "rules"), ...) {
  # Capture all arguments
  args <- list(...)
  verbose <- if ("verbose" %in% names(args)) args$verbose else TRUE

  # If no arguments, nothing to validate
  if (length(args) == 0) {
    return(return(args))
  }

  # Validate each argument in order
  for (arg_name in names(args)) {
    # Skip unnamed arguments
    if (is.null(arg_name) || arg_name == "") {
      next
    }

    # Check if we have a rule for this argument (FOR DEBUGGIN ONLY)
    if (!arg_name %in% names(rules.mapping)) {
      if (verbose) {
        .warning.message(sprintf("No validation rule defined for '%s'", arg_name))
      }
      next
    }

    # Build the validation rule
    # WARNING: arguments must be passed in order of priority
    value <- args[[arg_name]]
    rule <- do.call(rules.mapping[[arg_name]], c(list(x = value), args))

    # Validate the argument (para no primeiro erro)
    # update in case of any change
    args[[arg_name]] <- do.call(.validate.arg, c(list(value = value, name = arg_name), rule))
  }

  # Update environment if requested
  if (!is.null(envir)) {
    list2env(args, envir = envir)
    return(invisible(TRUE))
  }

  # Otherwise return the validated arguments
  return(args)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Check and format predictor matrix X #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.X <- function(X, n_obs = NULL, n_feat = NULL) {
  .rule.X(n_obs = n_obs, n_feat = n_feat)$custom(x = X)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Check and format the response y #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.check.y <- function(y) {
  .rule.y()$custom(x = y)
}
