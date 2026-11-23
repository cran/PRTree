## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 4.2,
  fig.align = "center"
)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

## ----load-pkg-----------------------------------------------------------------
library(PRTree)

## ----basic-fit----------------------------------------------------------------
set.seed(1234)
X <- matrix(runif(200, 0, 10), ncol = 1)
y <- cos(X[, 1]) + rnorm(200, 0, 0.08)

# Fit a PRTree using all observations for training (prop_hold = 0)
fit <- pr_tree(y, X, max_terminal_nodes = 9, prop_hold = 0)
summary(fit)

## ----plot-1d-fit--------------------------------------------------------------
ord <- order(X[, 1])
plot(
  X[ord, 1], y[ord],
  pch = 16, cex = 0.6, col = "gray50",
  xlab = "x", ylab = "y",
  main = "PRTree Smooth Fit to cos(x)"
)
lines(X[ord, 1], fit$yhat[ord], col = "steelblue", lwd = 2)
lines(X[ord, 1], cos(X[ord, 1]), col = "firebrick", lwd = 1.5, lty = 2)
legend(
  "topright",
  legend = c("Observed", "PRTree fit", "True cos(x)"),
  col = c("gray50", "steelblue", "firebrick"),
  pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1.5),
  bty = "n"
)

## ----plot-tree, fig.height = 4.6----------------------------------------------
plot_tree(fit, max_bins = 10)

## ----diag-plots, eval = has_ggplot2, fig.height = 5---------------------------
plot(fit, y, which = c(1, 2, 4, 5), ncol = 2)

## ----predict-basic------------------------------------------------------------
X_new <- matrix(c(1.0, 3.14, 6.28, 9.0), ncol = 1)
pred <- predict(fit, newdata = X_new, complete = TRUE)

round(pred$yhat, 4)
round(rowSums(pred$P), 6)

## ----sigma-and-dist-----------------------------------------------------------
set.seed(42)
X2 <- matrix(runif(200 * 2, 0, 10), ncol = 2)
colnames(X2) <- c("x1", "x2")
y2 <- sin(X2[, 1]) + 0.4 * X2[, 2] + rnorm(200, 0, 0.15)

# Inspect an automatically scaled sigma grid
grid <- expand_sigma_grid(X2, grid_size = 5, tiny_sigma = 1e-20, verbose = FALSE)
grid

# Fit using Student's t distribution (df = 4) and the custom grid
ctrl <- pr_tree_control(
  dist = "t",
  df = 4,
  sigma_grid = grid,
  max_terminal_nodes = 8,
  prop_hold = 0.2
)
fit_t <- pr_tree(y2, X2, control = ctrl)
fit_t

## ----fixed-sigma--------------------------------------------------------------
ctrl_crisp <- pr_tree_control(sigma_grid = 0, prop_hold = 0, max_terminal_nodes = 6)
fit_crisp <- pr_tree(y2, X2, control = ctrl_crisp)
fit_crisp$sigma

## ----missing-data-------------------------------------------------------------
set.seed(2026)
X_miss <- X2
X_miss[sample(nrow(X_miss), 25), 1] <- NA
X_miss[sample(nrow(X_miss), 20), 2] <- NA

# Compare the three missing-data probability strategies
mse_by_fill <- sapply(0:2, function(ft) {
  m <- pr_tree(
    y2, X_miss,
    fill_type = ft,
    proxy_crit = "both",
    max_terminal_nodes = 8,
    prop_hold = 0.2,
    grid_size = 4
  )
  m$mse["validation"]
})
names(mse_by_fill) <- paste0("fill_type=", 0:2)
round(mse_by_fill, 4)

## ----predict-missing----------------------------------------------------------
fit_na <- pr_tree(
  y2, X_miss,
  fill_type = 2,
  max_terminal_nodes = 8,
  prop_hold = 0,
  grid_size = 4
)

X_new_na <- rbind(
  c(2.0, 5.0),
  c(NA, 5.0),
  c(2.0, NA),
  c(NA, NA)
)
colnames(X_new_na) <- colnames(X2)

predict(fit_na, newdata = X_new_na)

## ----split--------------------------------------------------------------------
set.seed(123)
is_na_row <- apply(is.na(X_miss), 1, any)
folds <- train_test_split(stratify = TRUE, is_NA = is_na_row, n_rep = 5)
summary(folds)

## ----cv-run-------------------------------------------------------------------
set.seed(123)
ctrl_cv <- pr_tree_control_cv(
  method = "kfold",
  n_rep = 5,
  prop_valid = 0.2,
  stratify = TRUE,
  max_terminal_nodes = 6,
  grid_size = 3
)

cv_res <- pr_tree_cv(y2, X_miss, control_cv = ctrl_cv, verbose = FALSE)
summary(cv_res)

## ----cv-plot, eval = has_ggplot2, fig.height = 5------------------------------
plot(cv_res)

