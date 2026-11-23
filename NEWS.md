# PRTree 1.1.1 (YYYY-MM-DD)

## Bug fixes

- Fixed a bug in the handling of `sigma_grid`: a scalar or a vector supplied through `pr_tree_control()` (e.g., `sigma_grid = 0`) was converted to a `1 x 1` matrix before the number of predictors was known. The Fortran routine then read beyond the supplied values, and `model$sigma` contained uninitialized memory. The grid is now expanded to `n_feat` columns once `X` is available.
- Grids with an incompatible number of columns (a vector with a length different from the number of predictors, or a matrix with the wrong number of columns) were silently accepted. They now stop with an informative error. Also fixed the malformed `sprintf()` call that produced the error message for vectors of the wrong length.
- Added the internal function `.check.sigma_grid()`, following the `.check.X()` and `.check.y()` pattern, and a consistency check of `sigma_grid` before calling the Fortran routine.
- Explicitly allocated `p_zero` and `idx_node` before assignment in `src/prtree.f90` and initialized `best_train_mse` in `src/Rcalls.f90` to resolve `-Wuninitialized` Fortran compiler warnings (including on macOS / M1mac).

## Documentation and maintenance

- Added the published DOI for the theoretical paper (`10.1016/j.apm.2026.117450`) to `DESCRIPTION`, `inst/CITATION`, and documentation, while keeping the R package paper reference pointing to `arXiv:2510.03634`.
- Added `README.md`, `inst/CITATION`, an introductory vignette (`vignette("PRTree-intro")`), and a `testthat` test suite, and moved `inst/news.md` to `NEWS.md` in the package root.
- Fixed minor documentation inconsistencies (typos, `prop_hold` and `n_rep` parameter names in examples/details, and restoring graphical parameters in `plot_tree()` and examples).

# PRTree 1.1.0 (2026-07-13)

## New features

- Added support for cross-validation through new functions for model fitting, evaluation, and result summaries (`pr_tree_cv()`, `pr_tree_control_cv()`, `train_test_split()`).

## Improvements

- Refactored the internal R and Fortran codebase to improve maintainability and extensibility while preserving backward compatibility.
- Improved the documentation of several functions and parameters.
- Corrected the documented default value of `n_min`.

## Changes

- Renamed `perc_test` to `prop_hold` to better reflect its role as the proportion of observations held out for validation or testing, depending on the `sigma_grid` setting.
- Renamed `perc_x` to `prop_x` for consistency with the package naming convention.

## Bug fixes

- Fixed an error in the Fortran implementation for non-Gaussian distributions. The standardized variable is now correctly computed as `(v - mu) / sigma` instead of `v * sigma + mu`.

# PRTree 1.0.3 (2026-02-18)

- Fixed an issue in `pr_tree_control()` where distribution parameters (e.g., `df` for t-distribution) were not properly recognized when creating a control list first and then passing it to `pr_tree()`. Direct assignment to `pr_tree()` via `...` worked correctly. Parameters now work in both scenarios.
- Added explicit type conversion for `y` and `X` to ensure compatibility with Fortran code.

# PRTree 1.0.2 (2025-10-26)

- Fixed a bug where the `n_train` argument wasn't being processed correctly for a particular scenario, causing `NULL` to be passed to the building function.
- Fixed a bug that caused an incorrect computation of the total number of missing values in Fortran, which led to crashes in specific scenarios.

# PRTree 1.0.1 (2025-10-18)

- Fixed a bug where the `fill_type` argument wasn't being processed correctly, causing `NULL` to be passed to the prediction function.

# PRTree 1.0.0 (2025-10-09)

Major update: complete architectural overhaul of the package, introducing significant new features, performance improvements, and enhanced stability.

## New features

- **Complete code restructuring** into modular Fortran components for improved maintainability and a cleaner separation of concerns.
- **Enhanced missing value handling** with three new strategies via the `fill_type` parameter:
  - `0`: Uniform probability.
  - `1`: Probability = 1 if all non-missing features are in the region.
  - `2`: Probability based only on non-missing features.
- **New proxy criteria** for assigning observations with missing values during a split, controlled by the `crit` parameter:
  - `mean`: Maximize the difference in means between child nodes.
  - `var`: Maximize the between-node variability.
  - `both`: Combine both criteria.
- **Controllable verbosity** for logging the tree-building process, allowing for different levels of detail from quiet to full debug output.
- **Enhanced parametrization**:
  - Each variable can now have a different `sigma`. The sigma values are passed via `sigma_grid` (replaces `sigmas` in the old version).
  - Added the `grid_size` parameter for sigma grid generation.
  - The algorithm now performs validation to find the best `sigma` value. The size of the validation sample is controlled via `perc_test` (later renamed to `prop_hold` in `1.1.0`).
- **Multiple split candidates** can now be evaluated via the new `n_candidates` parameter.
- **Two search modes** for optimal splits:
  - Node-local candidate search (`by_node = TRUE`).
  - Global candidate search (`by_node = FALSE`).
- **Output enhancements**:
  - Added `XRegion` output showing the final region assignments for each observation.
  - Improved `NA` tracking in the output.

## Algorithm improvements

- **Numerical stability**: Switched to LAPACK's `DGELSD` for solving least squares problems. This routine uses a Singular Value Decomposition (SVD) method, providing robust minimum-norm solutions even for rank-deficient matrices.
- **Performance**: Significantly more efficient due to:
  - A more optimized grid search implementation.
  - Reduced data copying and better memory management with controlled allocations/deallocations.
  - An optimized Fortran interface.
  - Efficient probability calculations that focus only on features with finite bounds in each region.
  - Incremental updates after splits to avoid full recomputations.
- **Code quality**: Better documentation, a more consistent coding style, and enhanced parameter validation.
- **Expanded functionality**: Supports more complex use cases, handles edge conditions more effectively, and provides more configurable behavior.

## Backward-incompatible changes

- The old `Iindep` parameter has been removed (it was unused in the previous version).
- The default value for `max_depth` is now `max_terminal_nodes - 1`.
- The prediction function now returns only `yhat` by default. The full probability matrix `P` can be obtained with `complete = TRUE`. The `newdata` field has been removed from the output list.
- Added the `XRegion` field to the output.
- Requires additional parameters for enhanced features.

# PRTree 0.1.3 (2025-05-22)

- **Bug fix**: Fixed a critical bug in the `predict` function where the probability matrix `P` was not computed correctly due to a type mismatch (`0` instead of `0L` for a parameter), which resulted in incorrect predictions.
- Last stable release based on the original `base.f90` architecture (basic probabilistic regression tree functionality, single split candidate evaluation, and simple uniform-probability missing value handling).

# PRTree 0.1.2 (2024-09-28)

- **Bug fix**: Fixed a bug introduced in version `0.1.1`.

# PRTree 0.1.1 (2024-09-14)

- Addressed a compilation issue with the `flang-19` compiler.

# PRTree 0.1.0 (2024-01-15)

- Initial release of the `PRTree` package.
- Provided basic probabilistic regression tree functionality based on the `base.f90` implementation, featuring single split candidate evaluation and simple missing value handling (uniform probabilities).
