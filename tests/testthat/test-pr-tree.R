test_that("pr_tree fits a basic model and returns consistent components", {
  set.seed(1234)
  X <- matrix(runif(150, 0, 10), ncol = 1)
  y <- cos(X[, 1]) + rnorm(150, 0, 0.05)

  fit <- pr_tree(y, X, max_terminal_nodes = 6, prop_hold = 0, grid_size = 4)

  expect_s3_class(fit, "prtree")
  expect_equal(fit$n_obs, 150L)
  expect_equal(fit$n_feat, 1L)
  expect_length(fit$yhat, 150L)
  expect_true(all(is.finite(fit$yhat)))
  expect_equal(nrow(fit$P), 150L)
  expect_equal(ncol(fit$P), length(fit$gamma))
  expect_equal(unname(rowSums(fit$P)), rep(1, 150L), tolerance = 1e-6)
  expect_true(all(fit$P >= -1e-12 & fit$P <= 1 + 1e-12))
  expect_true(fit$mse["train"] < 0.1)
})

test_that("predict.prtree returns expected vector and complete output", {
  set.seed(1234)
  X <- matrix(runif(120 * 2, 0, 10), ncol = 2)
  y <- sin(X[, 1]) + 0.5 * X[, 2] + rnorm(120, 0, 0.1)

  fit <- pr_tree(y, X, max_terminal_nodes = 5, prop_hold = 0, grid_size = 3)
  X_new <- matrix(runif(20 * 2, 0, 10), ncol = 2)

  pred_vec <- predict(fit, newdata = X_new)
  expect_true(is.numeric(pred_vec))
  expect_length(pred_vec, 20L)
  expect_true(all(is.finite(pred_vec)))

  pred_full <- predict(fit, newdata = X_new, complete = TRUE)
  expect_type(pred_full, "list")
  expect_named(pred_full, c("yhat", "P"))
  expect_equal(pred_full$yhat, pred_vec)
  expect_equal(nrow(pred_full$P), 20L)
  expect_equal(ncol(pred_full$P), length(fit$gamma))
  expect_equal(rowSums(pred_full$P), rep(1, 20L), tolerance = 1e-6)

  # Predicting on training data matches fitted values
  pred_train <- predict(fit, newdata = X)
  expect_equal(as.numeric(pred_train), as.numeric(fit$yhat), tolerance = 1e-6)
})

test_that("single-terminal-node tree works in fitting and prediction", {
  set.seed(99)
  X <- matrix(runif(60 * 2, 0, 5), ncol = 2)
  y <- rnorm(60, mean = 3, sd = 0.2)

  fit <- pr_tree(y, X, max_terminal_nodes = 1, prop_hold = 0, sigma_grid = 0.5)
  expect_length(fit$gamma, 1L)
  expect_equal(ncol(fit$P), 1L)

  pred <- predict(fit, newdata = X[1:5, , drop = FALSE], complete = TRUE)
  expect_equal(pred$yhat, rep(fit$gamma, 5L))
  expect_equal(dim(pred$P), c(5L, 1L))
  expect_equal(as.numeric(pred$P), rep(1, 5L))
})

test_that("non-Gaussian distributions work via pr_tree_control and ...", {
  set.seed(2026)
  X <- matrix(runif(100 * 2, 0, 5), ncol = 2)
  y <- sin(X[, 1]) + 0.3 * X[, 2] + rnorm(100, 0, 0.1)

  # Student-t: control list vs direct ... arguments (regression for v1.0.3)
  ctrl_t <- pr_tree_control(
    dist = "t", df = 4, sigma_grid = c(0.5, 0.5),
    prop_hold = 0, max_terminal_nodes = 4
  )
  fit_t1 <- pr_tree(y, X, control = ctrl_t)
  fit_t2 <- pr_tree(
    y, X, dist = "t", df = 4, sigma_grid = c(0.5, 0.5),
    prop_hold = 0, max_terminal_nodes = 4
  )
  expect_equal(fit_t1$dist$dist_name, "t")
  expect_equal(fit_t1$dist$df, 4)
  expect_equal(fit_t1$yhat, fit_t2$yhat, tolerance = 1e-10)
  expect_equal(unname(rowSums(fit_t1$P)), rep(1, 100L), tolerance = 1e-6)

  # Log-normal and Gamma
  fit_ln <- pr_tree(
    y, X, dist = "lnorm", sdlog = 1, sigma_grid = c(0.5, 0.5),
    prop_hold = 0, max_terminal_nodes = 4
  )
  expect_equal(fit_ln$dist$dist_name, "lnorm")
  expect_true(all(is.finite(fit_ln$yhat)))
  expect_equal(unname(rowSums(fit_ln$P)), rep(1, 100L), tolerance = 1e-6)

  fit_ga <- pr_tree(
    y, X, dist = "gamma", shape = 2, sigma_grid = c(0.5, 0.5),
    prop_hold = 0, max_terminal_nodes = 4
  )
  expect_equal(fit_ga$dist$dist_name, "gamma")
  expect_true(all(is.finite(fit_ga$yhat)))
  expect_equal(unname(rowSums(fit_ga$P)), rep(1, 100L), tolerance = 1e-6)
})

test_that("by_node = TRUE and by_node = FALSE both run cleanly", {
  set.seed(7)
  X <- matrix(runif(100 * 2, 0, 5), ncol = 2)
  y <- X[, 1]^2 - 2 * X[, 2] + rnorm(100, 0, 0.2)

  fit_global <- pr_tree(y, X, by_node = FALSE, max_terminal_nodes = 5, prop_hold = 0, grid_size = 3)
  fit_local <- pr_tree(y, X, by_node = TRUE, max_terminal_nodes = 5, prop_hold = 0, grid_size = 3)

  expect_s3_class(fit_global, "prtree")
  expect_s3_class(fit_local, "prtree")
  expect_true(all(is.finite(fit_global$yhat)))
  expect_true(all(is.finite(fit_local$yhat)))
})
