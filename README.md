# PRTree: Probabilistic Regression Trees

[![CRAN status](https://www.r-pkg.org/badges/version/PRTree)](https://CRAN.R-project.org/package=PRTree)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI (Theoretical Paper)](https://img.shields.io/badge/DOI-10.1016%2Fj.apm.2026.117450-blue)](https://doi.org/10.1016/j.apm.2026.117450)
[![arXiv (Package Paper)](https://img.shields.io/badge/arXiv-2510.03634-b31b1b.svg)](https://doi.org/10.48550/arXiv.2510.03634)

## Overview

**PRTree** provides an efficient R implementation of **Probabilistic Regression Trees (PRTrees)** with built-in handling of missing predictor values. Unlike standard regression trees that rely on hard binary splits, PRTrees assign observations to terminal regions smoothly via probability distributions, yielding continuous regression surfaces.

The core tree-building and prediction routines are implemented in **Fortran** (using LAPACK's `DGELSD` for numerically stable least-squares estimation), offering high computational efficiency alongside a clean S3 interface in R.

## Key Features

- **Smooth probabilistic splits:** Supports Gaussian (`"norm"`), log-normal (`"lnorm"`), Student's *t* (`"t"`), and gamma (`"gamma"`) distributions for computing region membership probabilities.
- **Native missing-data handling:** Accommodates `NA` values directly during both tree construction and prediction without requiring prior imputation, with three probability-filling strategies (`fill_type`) and configurable proxy split criteria (`proxy_crit`).
- **Automatic or custom smoothing grids:** Supports per-feature smoothing parameters (`sigma_grid`) generated automatically from training data scales or customized via `expand_sigma_grid()`.
- **Cross-validation tools:** Built-in *k*-fold and Monte Carlo cross-validation (`pr_tree_cv()`, `pr_tree_control_cv()`, and `train_test_split()`) with optional missing-value stratification.
- **Diagnostics and visualization:** Standard S3 methods (`print()`, `summary()`, `predict()`, `plot()`) plus tree visualization with terminal-node probability heatmaps (`plot_tree()`).

## Installation

Install the released version from CRAN:

```r
install.packages("PRTree")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("tsprass/PRTree")
```

## Quick Example

```r
library(PRTree)

# Simulate noisy data
set.seed(1234)
X <- matrix(runif(200, 0, 10), ncol = 1)
y <- cos(X) + rnorm(200, 0, 0.05)

# Fit a Probabilistic Regression Tree
fit <- pr_tree(y, X, max_terminal_nodes = 9, prop_hold = 0)
summary(fit)

# Predict on new observations (including missing values)
X_new <- matrix(c(2.5, 5.0, NA), ncol = 1)
predict(fit, newdata = X_new)

# Diagnostic plots and tree visualization
plot(fit, y, which = c(1, 2, 4, 5), ncol = 2)
plot_tree(fit, max_bins = 10)
```

Cross-validation example:

```r
set.seed(123)
X <- matrix(runif(1000, 0, 10), ncol = 2)
y <- 2 * sin(X[, 1]) + 0.5 * X[, 2] + rnorm(500, 0, 0.2)

ctrl_cv <- pr_tree_control_cv(
  method = "montecarlo",
  n_rep = 10,
  prop_test = 0.2,
  prop_valid = 0.2
)
cv_res <- pr_tree_cv(y, X, control_cv = ctrl_cv)
summary(cv_res)
plot(cv_res)
```

## Main Functions

| Function | Description |
| :--- | :--- |
| `pr_tree()` | Fit a Probabilistic Regression Tree model |
| `predict()` | Generate predictions (and optionally probability matrices) for new data |
| `pr_tree_control()` | Configure tree-building, distribution, and missing-data parameters |
| `expand_sigma_grid()` | Build custom smoothing parameter (`sigma`) grids |
| `pr_tree_cv()` | Perform *k*-fold or Monte Carlo cross-validation |
| `pr_tree_control_cv()` | Configure cross-validation settings |
| `train_test_split()` | Create random or missing-value stratified splits and folds |
| `plot_tree()` | Visualize the fitted tree structure and terminal probability distributions |

## References

- **Theoretical framework (missing data in PRTrees):**
  Prass, T. S., Neimaier, A. S., & Pumi, G. (2026). *Handling Missing Data in Probabilistic Regression Trees*. Applied Mathematical Modelling, 117450. <https://doi.org/10.1016/j.apm.2026.117450>

- **Software and implementation paper:**
  Prass, T. S., Neimaier, A. S., & Pumi, G. (2025). *PRTree: An R Package for Probabilistic Regression Trees with Built-in Missing Data Handling*. arXiv:2510.03634. <https://doi.org/10.48550/arXiv.2510.03634>

- **Original PRTree methodology:**
  Alkhoury, S., Devijver, E., Clausel, M., Tami, M., Gaussier, E., & Oppenheim, G. (2020). *Smooth and Consistent Probabilistic Regression Trees*. Advances in Neural Information Processing Systems (NeurIPS), 33, 11345–11355.

## Citation

To cite `PRTree` in publications, run `citation("PRTree")` in R.

## License

GPL (>= 3)
