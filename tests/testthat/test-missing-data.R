test_that("missing values are handled across all fill_type strategies", {
  set.seed(101)
  n <- 120
  X <- matrix(runif(n * 3, 0, 10), ncol = 3)
  y <- sin(X[, 1]) + 0.5 * X[, 2] - 0.3 * X[, 3] + rnorm(n, 0, 0.1)

  # Inject missing values in predictors
  na_rows <- c(1, 5, 10, 20, 35, 50, 75, 100)
  X[na_rows[1:4], 1] <- NA
  X[na_rows[5:6], 2] <- NA
  X[na_rows[7:8], 3] <- NA

  X_new <- rbind(
    c(2.0, 5.0, 3.0),
    c(NA, 5.0, 3.0),
    c(2.0, NA, NA),
    c(NA, NA, NA)
  )

  for (ft in c(0L, 1L, 2L)) {
    fit <- pr_tree(
      y, X,
      fill_type = ft,
      max_terminal_nodes = 5,
      prop_hold = 0,
      sigma_grid = c(0.5, 0.5, 0.5)
    )

    expect_equal(unname(fit$fill_type), ft)
    expect_true(all(is.finite(fit$yhat)))
    expect_equal(unname(rowSums(fit$P)), rep(1, n), tolerance = 1e-6)
    expect_true(all(is.na(fit$XRegion[na_rows, "TRUE"])))
    expect_true(all(!is.na(fit$XRegion[, "Internal"])))

    pred <- predict(fit, newdata = X_new, complete = TRUE)
    expect_true(all(is.finite(pred$yhat)))
    expect_equal(rowSums(pred$P), rep(1, nrow(X_new)), tolerance = 1e-6)
    expect_true(all(pred$P >= -1e-12 & pred$P <= 1 + 1e-12))
  }
})

test_that("all proxy_crit options work with missing data", {
  set.seed(202)
  n <- 100
  X <- matrix(runif(n * 2, 0, 10), ncol = 2)
  y <- 2 * sin(X[, 1]) + X[, 2] + rnorm(n, 0, 0.2)
  X[sample(n, 15), 1] <- NA
  X[sample(n, 10), 2] <- NA

  for (crit in c("mean", "var", "both")) {
    fit <- pr_tree(
      y, X,
      proxy_crit = crit,
      max_terminal_nodes = 5,
      prop_hold = 0,
      sigma_grid = c(0.4, 0.4)
    )
    expect_s3_class(fit, "prtree")
    expect_true(all(is.finite(fit$yhat)))
    expect_equal(unname(rowSums(fit$P)), rep(1, n), tolerance = 1e-6)
  }
})
