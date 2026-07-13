# Last revision: July/2026

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Collects all arguments from the parent function call, #----
# including default values for arguments that were not supplied.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.get.input <- function(...) {
  # Capture explicit arguments from parent frame
  parent <- parent.frame()
  var_names <- ls(envir = parent)

  # Get explicit arguments
  control <- mget(var_names, envir = parent)

  # Add ... arguments
  dots <- list(...)
  if (length(dots) > 0) {
    control <- c(control, dots)
  }

  return(control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: set storage mode #----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.set.storage_mode <- function(x, mode) {
  if (!is.null(x)) {
    storage.mode(x) <- mode
  }
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal. Process the control list #----
# 1) Update the tree building parameters (and validate)
#    Overrides defaults with user-provided values
#    Adds extra information needed to further processing
#    Returns with updated = TRUE, so updating will be skiped next time
# 2) Updates the cv parameters
#    Adds extra information needed to further processing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @importFrom utils modifyList
.update.control <- function(y, X, control, params = list(), is.cv = FALSE,
                            need_idx = !is.cv, need_sigma = !is.cv) {
  # make sure control is a list
  if (!is.list(control)) .stop.message("  `control` must be a list")

  # Add basic info from data
  control$n_obs <- length(y)
  control$n_feat <- ncol(X)

  # Merge with params (params take precedence)
  all_params <- modifyList(control, params)

  # Validate everything at once
  validated <- .pr_tree_control(all_params, is.cv = is.cv)
  control[names(validated)] <- validated

  # get the indexes for the training sample.
  # If indexes are provided, only sort them, otherwhise create a sample using
  # the stratified sampling procedure (if prop_test > 0)
  # skip this step in CV and compute inside the loop
  if (need_idx) {
    train <- .sample.idx(
      idx_train = control$idx_train,
      is_NA = apply(is.na(X), 1, any),
      n_obs = control$n_obs,
      prop_test = control$prop_hold,
      stratify = TRUE
    )
    control[names(train)] <- train
  }

  # Generate sigma_grid if necessary, using a grid based on the
  # variances of X (validation set is not used to create the grid).
  #
  # For CV:
  #  - skip this step and compute inside the loop when grid_size > 0
  #  - build sigma_grid when grid_size = 0
  if (need_sigma && is.null(control$sigma_grid)) {
    control$sigma_grid <- .expand_sigma_grid(
      X = X, n_feat = ncol(X),
      idx_train = control$idx_train,
      grid_size = control$grid_size,
      min_mult = control$min_mult,
      max_mult = control$max_mult,
      tiny_sigma = control$tiny_sigma,
      verbose = control$verbose
    )
  } else if (is.cv && is.null(control$sigma_grid) &&
            !is.null(control$tiny_sigma) && control$grid_size < 1) {
    control$sigma_grid <- matrix(control$tiny_sigma, ncol = control$n_feat)
  }

  # Final grid size
  control$grid_size_final <- .grid_size_final(
    sigma_grid = control$sigma_grid,
    grid_size = control$grid_size,
    tiny_sigma = control$tiny_sigma
  )

  # Debug default
  if (is.null(control$debug)) {
    control$debug <- .pr.meta_list$debug$default
  }

  # CV-specific computed fields
  if (is.cv) {
    # null_sigma: TRUE if sigma_grid = NULL
    #  in which case the sigma_grid will be computed during
    #  the validation process
    control$null_sigma <- is.null(control$sigma_grid)

    # multiple_sigma: TRUE if grid_size > 1 (sigma_grid can be know or NULL)
    # - if TRUE, sigma_grid need to be updated in step 2 of fit.tree
    # - if FALSE, only one value of sigma is provided and prop_valid = 0
    #   by default (.rule.prop.valid)
    control$multiple_sigma <- control$grid_size_final > 1

    # track_sigma: TRUE if sigma is not a fixed vector
    # (will be used to print/save sigma related results)
    control$track_sigma <- control$null_sigma || control$multiple_sigma

    # simplified: TRUE when there are multiple sigmas and prop_valid = 0
    # (will be used to colect the mse value from the output)
    null_valid <- is.null(control$prop_valid)
    if (null_valid) {
      control$simplified <- FALSE
    } else {
      control$simplified <- control$track_sigma && !(control$prop_valid > 0)
    }

    # track_test: TRUE if a test set exists
    # will be used to print/update results related to testing
    if (is.null(control$prop_test)) {
      control$track_test <- FALSE
    } else {
      control$track_test <- control$prop_test > 0
    }

    # update_final: no reason to update if prop_valid = 0
    if(control$simplified) control$update_final = FALSE

		# output_mode: CV/light output
    need_output <- !control$update_final && control$track_test
    control$output_mode <- if(need_output) 1L else 0L
  } else {
		control$output_mode <- 1L
	}

  return(control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: Get the final grid size #----
# based on the provided sigma_grid and tiny_sigma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.grid_size_final <- function(sigma_grid, grid_size, tiny_sigma) {
  if (!is.null(sigma_grid)) {
    return(nrow(sigma_grid))
  }
  if (is.null(tiny_sigma)) {
    return(grid_size)
  }
  return(grid_size + 1)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: update idx_train in cv #----
# Update idx_train and n_train consistently during cv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.update.idx_train <- function(control, new_idx, n_obs) {
  storage.mode(new_idx) <- "integer"
  control$idx_train <- new_idx
  control$n_train <- length(new_idx)
  control$n_obs <- as.integer(n_obs)
  return(control)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal: update sigma in cv #----
# Update sigma_grid and grid_size consistently durign cv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.update.sigma_grid <- function(control, sigma_grid) {
  control$sigma_grid <- sigma_grid
  control$grid_size_final <- nrow(sigma_grid)
  return(control)
}
