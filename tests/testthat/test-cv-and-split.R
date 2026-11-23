test_that("train_test_split creates valid simple, stratified, and k-fold splits", {
  set.seed(123)

  # Simple split
  sp1 <- train_test_split(n_obs = 100, prop_test = 0.3)
  expect_s3_class(sp1, "idx.split")
  expect_equal(sp1$n_train, 70L)
  expect_equal(sp1$n_test, 30L)
  expect_length(intersect(sp1$idx_train, sp1$idx_test), 0L)
  expect_equal(sort(c(sp1$idx_train, sp1$idx_test)), 1:100)

  # Stratified split by missing values
  is_na <- rep(c(TRUE, FALSE), times = c(20, 80))
  sp2 <- train_test_split(prop_test = 0.2, stratify = TRUE, is_NA = is_na)
  expect_true(sp2$stratified)
  expect_equal(sp2$missing_stats$total_na, 20L)
  expect_equal(sp2$missing_stats$train_na + sp2$missing_stats$test_na, 20L)
  expect_equal(sp2$missing_stats$test_na, 4L)

  # Stratified k-fold partition
  folds <- train_test_split(stratify = TRUE, is_NA = is_na, n_rep = 5)
  expect_s3_class(folds, "idx.split")
  expect_equal(folds$n_rep, 5L)
  expect_length(folds$fold_idx, 100L)
  expect_equal(sort(unique(folds$fold_idx)), 1:5)
  expect_equal(sum(folds$fold_stats$n_na), 20L)
})

test_that("pr_tree_cv runs Monte Carlo, k-fold, only_sigma, and fixed sigma", {
  set.seed(321)
  n <- 100
  X <- matrix(runif(n * 2, 0, 5), ncol = 2)
  y <- sin(X[, 1]) + 0.4 * X[, 2] + rnorm(n, 0, 0.15)

  # Monte Carlo CV
  ctrl_mc <- pr_tree_control_cv(
    method = "montecarlo", n_rep = 3, prop_test = 0.2, prop_valid = 0.2,
    max_terminal_nodes = 4, grid_size = 2
  )
  cv_mc <- pr_tree_cv(y, X, control_cv = ctrl_mc, verbose = FALSE)
  expect_s3_class(cv_mc, "prtree.cv")
  expect_equal(dim(cv_mc$sigma_matrix), c(3L, 2L))
  expect_length(cv_mc$rmse_by_rep$validation, 3L)
  expect_length(cv_mc$rmse_by_rep$test, 3L)
  expect_true(all(is.finite(cv_mc$rmse_by_rep$test)))

  # k-fold CV
  ctrl_kf <- pr_tree_control_cv(
    method = "kfold", n_rep = 3, prop_valid = 0.2,
    max_terminal_nodes = 4, grid_size = 2
  )
  cv_kf <- pr_tree_cv(y, X, control_cv = ctrl_kf, verbose = FALSE)
  expect_s3_class(cv_kf, "prtree.cv")
  expect_equal(dim(cv_kf$sigma_matrix), c(3L, 2L))
  expect_length(cv_kf$rmse_by_rep$test, 3L)

  # Sigma selection only (no test set)
  ctrl_os <- pr_tree_control_cv(
    method = "kfold", n_rep = 3, only_sigma = TRUE,
    max_terminal_nodes = 4, grid_size = 2, verbose = FALSE
  )
  cv_os <- pr_tree_cv(y, X, control_cv = ctrl_os, verbose = FALSE)
  expect_null(cv_os$rmse_by_rep$test)
  expect_length(cv_os$rmse_by_rep$validation, 3L)

  # Fixed sigma (scalar expanded via pr_tree_control_cv)
  ctrl_fix <- pr_tree_control_cv(
    method = "kfold", n_rep = 3, sigma_grid = 0.5, prop_valid = 0,
    max_terminal_nodes = 4
  )
  cv_fix <- pr_tree_cv(y, X, control_cv = ctrl_fix, verbose = FALSE)
  expect_false(cv_fix$track_sigma)
  expect_length(cv_fix$rmse_by_rep$test, 3L)
})

test_that("print and summary S3 methods execute cleanly", {
  set.seed(11)
  X <- matrix(runif(80 * 2, 0, 5), ncol = 2)
  y <- cos(X[, 1]) + rnorm(80, 0, 0.1)

  ctrl <- pr_tree_control(max_terminal_nodes = 4, prop_hold = 0, grid_size = 2)
  ctrl_cv <- pr_tree_control_cv(n_rep = 2, max_terminal_nodes = 4, grid_size = 2)
  fit <- pr_tree(y, X, control = ctrl)
  cv <- pr_tree_cv(y, X, control_cv = ctrl_cv, verbose = FALSE)
  sp <- train_test_split(n_obs = 80, prop_test = 0.25)
  grid <- expand_sigma_grid(X, grid_size = 3, verbose = FALSE)

  expect_output(print(ctrl))
  expect_output(print(ctrl_cv))
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_output(print(cv))
  expect_output(print(summary(cv)))
  expect_output(print(sp))
  expect_output(print(summary(sp)))
  expect_output(print(grid))
})
