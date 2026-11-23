test_that("scalar sigma_grid via pr_tree_control expands to n_feat columns", {
  set.seed(123)
  X <- matrix(runif(120 * 3, 0, 10), ncol = 3)
  y <- sin(X[, 1]) + 0.5 * X[, 2] - 0.2 * X[, 3] + rnorm(120, 0, 0.1)

  ctrl_zero <- pr_tree_control(sigma_grid = 0, prop_hold = 0, max_terminal_nodes = 5)
  fit_zero <- pr_tree(y, X, control = ctrl_zero)
  expect_length(fit_zero$sigma, 3L)
  expect_equal(as.numeric(fit_zero$sigma), c(0, 0, 0))

  ctrl_val <- pr_tree_control(sigma_grid = 0.5, prop_hold = 0, max_terminal_nodes = 5)
  fit_val <- pr_tree(y, X, control = ctrl_val)
  expect_length(fit_val$sigma, 3L)
  expect_equal(as.numeric(fit_val$sigma), c(0.5, 0.5, 0.5))

  # Passing scalar directly through ... should give identical sigma
  fit_dots <- pr_tree(y, X, sigma_grid = 0.5, prop_hold = 0, max_terminal_nodes = 5)
  expect_equal(as.numeric(fit_dots$sigma), c(0.5, 0.5, 0.5))
})

test_that("vector and matrix sigma_grid with valid dimensions work properly", {
  set.seed(123)
  X <- matrix(runif(100 * 3, 0, 10), ncol = 3)
  y <- sin(X[, 1]) + 0.3 * X[, 2] + rnorm(100, 0, 0.1)

  # Vector of length n_feat
  vec_sigma <- c(0.2, 0.4, 0.6)
  ctrl_vec <- pr_tree_control(sigma_grid = vec_sigma, prop_hold = 0, max_terminal_nodes = 4)
  fit_vec <- pr_tree(y, X, control = ctrl_vec)
  expect_equal(as.numeric(fit_vec$sigma), vec_sigma)

  # Matrix with multiple candidate rows and n_feat columns
  mat_sigma <- rbind(c(0.1, 0.1, 0.1), c(0.5, 0.5, 0.5), c(1.0, 1.0, 1.0))
  fit_mat <- pr_tree(y, X, sigma_grid = mat_sigma, prop_hold = 0.2, max_terminal_nodes = 4)
  expect_length(fit_mat$sigma, 3L)
  expect_true(any(apply(mat_sigma, 1, function(r) isTRUE(all.equal(r, as.numeric(fit_mat$sigma))))))
})

test_that("incompatible sigma_grid dimensions throw informative errors", {
  set.seed(123)
  X <- matrix(runif(80 * 3, 0, 10), ncol = 3)
  y <- sin(X[, 1]) + rnorm(80, 0, 0.1)

  # Vector of wrong length (passed directly and via pr_tree_control)
  expect_error(
    pr_tree(y, X, sigma_grid = c(0.1, 0.2), prop_hold = 0),
    "columns"
  )
  ctrl_bad_vec <- pr_tree_control(sigma_grid = c(0.1, 0.2), prop_hold = 0)
  expect_error(
    pr_tree(y, X, control = ctrl_bad_vec),
    "columns"
  )

  # Matrix with wrong number of columns
  bad_mat <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)
  expect_error(
    pr_tree(y, X, sigma_grid = bad_mat, prop_hold = 0),
    "columns"
  )
  ctrl_bad_mat <- pr_tree_control(sigma_grid = bad_mat, prop_hold = 0)
  expect_error(
    pr_tree(y, X, control = ctrl_bad_mat),
    "columns"
  )
})

test_that("expand_sigma_grid builds expected dimensions and attributes", {
  set.seed(42)
  X <- matrix(runif(100 * 3, 0, 10), ncol = 3)
  colnames(X) <- c("A", "B", "C")

  g1 <- expand_sigma_grid(X, grid_size = 5, verbose = FALSE)
  expect_s3_class(g1, "sigma_grid")
  expect_equal(dim(g1), c(5L, 3L))
  expect_equal(colnames(g1), c("A", "B", "C"))
  expect_false(attr(g1, "has_tiny"))

  g2 <- expand_sigma_grid(n_feat = 3, grid_size = 4, tiny_sigma = 1e-20, verbose = FALSE)
  expect_equal(dim(g2), c(5L, 3L))
  expect_true(attr(g2, "has_tiny"))
  expect_equal(as.numeric(g2[1, ]), rep(1e-20, 3))

  g3 <- expand_sigma_grid(n_feat = 2, grid_size = 0, tiny_sigma = 1e-20, verbose = FALSE)
  expect_equal(dim(g3), c(1L, 2L))
  expect_equal(as.numeric(g3[1, ]), rep(1e-20, 2))
})
