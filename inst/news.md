# PRTree News

## Version 1.0.0 (Major Update)

This version is a complete architectural overhaul of the package, introducing significant new features, performance improvements, and enhanced stability.

### New Features

-   **Complete code restructuring** into modular Fortran components for improved maintainability and a cleaner separation of concerns.

-   **Enhanced missing value handling** with three new strategies via the `fill_type` parameter:

    -   `0`: Uniform probability.
    -   `1`: Probability = 1 if all non-missing features are in the region.
    -   `2`: Probability based only on non-missing features.

-   **New proxy criteria** for assigning observations with missing values during a split, controlled by the `crit` parameter:

    -   `mean`: Maximize the difference in means between child nodes.
    -   `var`: Maximize the between-node variability.
    -   `both`: Combine both criteria.

-   **Controllable verbosity** for logging the tree-building process, allowing for different levels of detail from quiet to full debug output.

-   **Enhanced Parametrization**:

    -   Each variable can now have a different `sigma`. The sigma values are passed via `sigma_grid` (replaces `sigmas` in the old version).
    -   Added the `grid_size` parameter for sigma grid generation.
    -   The algorithm now performs cross-validation to find the best `sigma` value. The size of the validation sample is controlled via `perc_test`.

-   **Multiple split candidates** can now be evaluated via the new `n_candidates` parameter.

-   **Two search modes** for optimal splits:

    -   Node-local candidate search (`by_node = TRUE`).
    -   Global candidate search (`by_node = FALSE`).

-   **Output Enhancements**:

    -   Added `XRegion` output showing the final region assignments for each observation.
    -   Improved NA tracking in the output.

### Algorithm Improvements

-   **Numerical Stability**: Switched to LAPACK's `DGELSD` for solving least squares problems. This routine uses a Singular Value Decomposition (SVD) method, providing robust minimum-norm solutions even for rank-deficient matrices.

-   **Performance**: The new version is significantly more efficient due to:

    -   A more optimized grid search implementation.
    -   Reduced data copying and better memory management with controlled allocations/deallocations.
    -   An optimized Fortran interface.
    -   Efficient probability calculations that focus only on features with finite bounds in each region.
    -   Incremental updates after splits to avoid full recomputations.

-   **Code Quality**: The code has been improved with better documentation, a more consistent coding style, and enhanced parameter validation.

-   **Expanded Functionality**: The new algorithm supports more complex use cases, handles edge conditions more effectively, and provides more configurable behavior.

### Backward-Incompatible Changes

-   The old `Iindep` parameter has been removed (it was unused in the previous version).
-   The default value for `max_depth` is now `max_terminal_node - 1`.
-   The prediction function now returns only `yhat` by default. The full probability matrix `P` can be obtained with `complete = TRUE`. The `newdata` field has been removed from the output list.
-   The new version adds `XRegion` field to the output
-   The new version requires additional parameters for enhanced features.

## Version 0.1.3 (Previous Stable)

-   Initial stable implementation (`base.f90`)
-   Basic probabilistic regression tree functionality
-   Single split candidate evaluation
-   Simple missing value handling (uniform probabilities)

# Version History

## PRTree 1.0.2 (Release Date: YYYY-MM-DD)

Fixed a bug where the `n_train` argument wasn't being processed correctly for a particular scenario, causing `NULL` to be passed to the building function.

Fixed a bug that caused an incorrect computation of the total number of missing values in Fortran, which led to crashes in specific scenarios. 

## PRTree 1.0.1 (Release Date: 2025-10-18)

Fixed a bug where the `fill_type` argument wasn't being processed correctly, causing `NULL` to be passed to the prediction function.


## PRTree 1.0.0 (2025-10-09)

-   Current Version

## PRTree 0.1.3 (2025-05-22)

-   **Bug Fix**: Fixed a critical bug in the `predict` function where the probability matrix `P` was not computed correctly due to a type mismatch (`0` instead of `0L` for a parameter), which resulted in incorrect predictions.

## PRTree 0.1.2 (2024-09-28)

-   **Bug Fix**: Fixed a bug introduced in version 0.1.1.

## PRTree 0.1.1 (2024-09-14)

-   Addressed a compilation issue with the `flang-19` compiler.

## PRTree 0.1.0 (2024-01-15)

-   Initial release of the `PRTree` package.

-   This version provided basic probabilistic regression tree functionality based on the `base.f90` implementation, featuring single split candidate evaluation and simple missing value handling (uniform probabilities).

-   Submission date: 2024-01-15
