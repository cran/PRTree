# convert the distribution name to the FORTRAN code
.get.dist <- function(dist, dist_pars) {
  pars <- dist_pars
  switch(dist,
    norm = {
      my_dist <- 1L
      par_dist <- 0
    },
    lnorm = {
      my_dist <- 2L
      par_dist <- pars$sdlog
    },
    t = {
      my_dist <- 3L
      par_dist <- pars$df
    },
    gamma = {
      my_dist <- 4L
      par_dist <- pars$shape
    }
  )
  if (is.null(par_dist)) {
    stop(paste0(
      "Required parameters for distribution '",
      dist,
      "' are missing with no default"
    ))
  }
  return(list(dist_name = dist, dist_code = my_dist, dist_pars = par_dist))
}


# extract the indexes for the testing set
.get.idx <- function(X, idx_train, perc_test) {
  # list to return the output
  out <- list()

  # if the index for the training sample is provided, ignore perc_test
  if (!is.null(idx_train)) {
    out$idx_train <- as.integer(sort(unique(idx_train)))
    n_train <- length(out$idx_train)
    return(out)
  }

  # If perc_test is provided, use perc_test to obtain the training sample
  if (perc_test > 0 && perc_test < 1) {
    # Stratified sampling to preserve missing value proportions, as per documentation
    has_NA <- apply(is.na(X), 1, any)

    # Sampling from the complete group
    idx_noNA <- which(!has_NA)
    n_noNA <- length(idx_noNA)
    n_train_noNA <- round(n_noNA * (1 - perc_test))
    if (n_train_noNA > 1) {
      train_idx_noNA <- sample(idx_noNA, n_train_noNA)
    } else {
      train_idx_noNA <- NULL
      n_train_noNA <- 0
    }

    # Sampling from the group with NA
    n_NA <- nrow(X) - n_noNA
    n_train_NA <- round(n_NA * (1 - perc_test))
    if (n_train_NA > 1) {
      idx_NA <- which(has_NA)
      train_idx_NA <- sample(idx_NA, n_train_NA)
    } else {
      train_idx_NA <- NULL
      n_train_NA <- 0
    }

    out$n_train <- as.integer(n_train_noNA + n_train_NA)
    out$idx_train <- as.integer(c(train_idx_noNA, train_idx_NA))
    return(out)
  }

  # Use all data for training if perc_test is 0
  n_obs <- nrow(X)
  out$n_train <- as.integer(n_obs)
  out$idx_train <- as.integer(1:n_obs)
  return(out)
}

.get.sigma.grid <- function(X, idx_train, sigma_grid, grid_size) {
  # if sigma_grid is not provided, compute it based on the input matrix X
  if (is.null(sigma_grid)) {
    # compute the sd by column
    ss <- apply(X[idx_train, , drop = FALSE], 2, function(x) sd(x, na.rm = TRUE))

    # grid_size
    if (is.null(grid_size)) grid_size <- 8

    # if grid_size is 1, uses the sd as estimates
    if (grid_size <= 1) {
      return(matrix(ss, nrow = 1))
    }

    # for grid_size > 1, creates a grid with multiples of the sd
    multipliers <- seq(from = 0, to = 2, length.out = grid_size + 1)[-1]
    sigma_grid <- outer(multipliers, ss, FUN = "*")
    return(sigma_grid)
  }

  # if sigma_grid is provided, check the compatibility
  n_feat <- ncol(X)
  if (!is.matrix(sigma_grid)) {
    lsi <- length(sigma_grid)
    if (lsi > 1 && lsi != n_feat) {
      stop(paste0(
        "The regressor matrix has ", n_feat, " columns\n",
        "Please, provide either 1 or ", n_feat, " values for sigma_grid", "\n"
      ))
    }
    sigma_grid <- matrix(sigma_grid, ncol = n_feat)
    return(sigma_grid)
  }

  lsi <- ncol(sigma_grid)
  if (lsi != n_feat) {
    stop(paste0(
      "The regressor matrix has ", n_feat, " columns\n",
      "The sigma_grid provided has ", lsi, " columns\n"
    ))
  }
  return(sigma_grid)
}
