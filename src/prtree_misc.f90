! Utility functions

module prtree_misc
  use prtree_types
  implicit none

contains
  !=======================================
  !   Utility functions
  !=======================================
  function sort(n, x) result(xout)
    !-------------------------------------------------------
    ! Sorts a real array and returns the sorted array.
    !
    ! ARGUMENTS
    !   n [in] : Integer, number of elements.
    !   x [in] : Real(dp), input array.
    !
    ! Added: July/2025
    !-------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), allocatable :: xout(:)
    integer :: i
    integer, allocatable :: idx(:)
    real(dp) :: tmp_val

    allocate(xout(n))
    ! make a copy to avoid changing the original data
    xout = x
    if (n < 2) return

    ! Use the more efficient Heapsort (revsortr) and reverse the result.
    allocate(idx(n))
    idx = [(i, i=1, n)]
    call revsortr(xout, idx, n)  ! Sorts descending

    ! Reverse in-place to get ascending order
    do i = 1, n / 2
      tmp_val = xout(i)
      xout(i) = xout(n - i + 1)
      xout(n - i + 1) = tmp_val
    end do
  end function sort

  subroutine sort_xy(n, x, y, idx)
    !-----------------------------------------------------------------------------------
    ! Sorts array 'x' in ascending order and permutes arrais 'y' and 'idx' accordingly.
    ! Uses revsort (a heapsort algorithm) and then reverses the result to get
    ! ascending order.
    !
    ! ARGUMENTS
    !   n   [in]    : Integer, the size of the arrays.
    !   x   [inout] : Real(dp) array to be sorted.
    !   y   [inout] : Real(dp) array to be permuted like x.
    !   idx [inout] : Integer, the permutted indexes
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: idx(n)
    real(dp), intent(inout) :: x(n), y(n)
    ! Local variables
    integer :: i
    integer, allocatable :: idx_n(:)
    real(dp) :: tmp_val
    integer :: tmp_idx

    if (n < 2) return

    ! 1. Create an index array [1, 2, ..., n]
    allocate(idx_n(n))
    idx_n = [(i, i=1, n)]

    ! 2. Sort x in descending order and permut idx to match.
    call revsortr(x, idx_n, n)

    ! 3. Reorder the original y and idx values according to the permuted index.
    y = y(idx_n)
    idx = idx(idx_n)

    ! 4. Reverse the arrays in-place to get the final ascending order.
    do i = 1, n / 2
       tmp_val = x(i)
       x(i) = x(n - i + 1)
       x(n - i + 1) = tmp_val
       tmp_val = y(i)
       y(i) = y(n - i + 1)
       y(n - i + 1) = tmp_val
       tmp_idx = idx(i)
       idx(i) = idx(n - i + 1)
       idx(n - i + 1) = tmp_idx
    end do
  end subroutine sort_xy

  function lsquare(tree, a) result(x)
    !-------------------------------------------------------------------------------------
    ! Use LAPACK's DGELSD to solve the least squares problem
    !   minimum || b - A * x||
    !
    ! ARGUMENTS
    !   tree[inout] : tree_model, contains data and workspaces
    !   A [in] : Real(dp), tree%n_train x (tree%n_tn + 1) matrix
    !
    ! RETURNS
    !   x : Real(dp), coefficients (size n vector)
    !
    ! Notes:
    ! - DGELSD computes the minimum-norm solution to a real linear least squares problem:
    !       Minimize 2-norm(| b - A*x |)
    !   using the singular value decomposition (SVD) of A with a divide-and-conquer method.
    !   It is suitable for matrices A that may be rank-deficient.
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: a_copy(:, :), b_copy(:, :)
    real(dp), allocatable :: s(:), work(:)
    integer, allocatable :: iwork(:)
    real(dp) :: rcond
    integer :: m, n, lda, ldb, info, rank, lwork
    real(dp) :: work_query(1)
    integer :: iwork_query(1)

    m = size(a, 1)   ! number of rows
    n = size(a, 2)   ! number of columns
    rcond = -1.0_dp  ! rcond < 0 means: use machine precision
    lda = m
    ldb = max(m, n)
    allocate(x(n))

    ! the acual daty does not matter for the workspace query,
    ! only the dimensions do. So we can allocate with any values.
    allocate(a_copy(lda, n), source = 0.0_dp)
    allocate(b_copy(ldb, 1), source = 0.0_dp)
    allocate(s(min(m, n)), source = 0.0_dp)

    ! Reallocate workspace only if n > lsq_n
    if (n > tree%lsq_n) then
      ! Workspace query
      lwork = -1
      info = 0
      work_query(1) = 0.0_dp
      iwork_query(1) = 1
      call dgelsd(m, n, 1, a_copy, lda, b_copy, ldb, s, rcond, rank, work_query, lwork, iwork_query, info)

      tree%lsq_n = n
      tree%lsq_lwork = max(1, int(work_query(1)))
      tree%lsq_liwork = max(1, iwork_query(1))
    end if

    allocate(work(tree%lsq_lwork))
    allocate(iwork(tree%lsq_liwork))

    ! Copy input data (DGELSD overwrites its input arrays)
    a_copy = a
    b_copy(1:m, 1) = tree%y_train
    if (ldb > m) b_copy(m + 1:ldb, 1) = 0.0_dp

    call dgelsd(m, n, 1, a_copy, lda, b_copy, ldb, s, rcond, rank, work, tree%lsq_lwork, iwork, info)

    ! The solution is returned in the first k rows of B
    if (info == 0) then
      x = b_copy(1:n,1)
    else
      x = 0.0_dp
      if (tree%printinfo >= log_base) then
        call labelpr("Warning: DGELSD failed in `lsquare`", -1)
        call intpr1("Info code:", -1, info)
      end if
    end if
  end function lsquare

  subroutine find_thresholds(min_obs, x, x_mask, f, thresholds)
    !------------------------------------------------------------------------------
    ! Finds candidate thresholds for splitting a node. A threshold is the midpoint
    ! between two unique, adjacent feature values. A split is only valid if it
    ! results in child nodes each having at least min_obs observations.
    !
    ! ARGUMENTS
    !   min_obs    [in]  : Integer, minimum observations per node.
    !   n_train    [in]  : Integer, number of observations.
    !   n_feat     [in]  : Integer, number of features.
    !   x          [in]  : Real(dp), feature values.
    !   x_mask     [in]  : Logical, mask indicating which observations to consider.
    !   f          [in]  : Integer, index of the feature for which to find thresholds.
    !   thresholds [out] : Type(info_thr), thresholds found for the feature.
    !
    ! Added July/2025
    ! Last revision: May, 2026
    !   Removed slicing and added x_mask to consider only the relevant observations.
    !   This allows us to avoid copying data and sorting the entire array when only ]
    !   a subset of observations is relevant for finding thresholds.
    !------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: min_obs
    real(dp), intent(in) :: x(:, :)
    logical, intent(in) :: x_mask(:)
    integer, intent(in) :: f
    type(info_thr), intent(out) :: thresholds
    real(dp), allocatable :: xs(:)
    integer :: i
    integer :: n_start, n_end, n_unique, n_consider

    ! Find the first and last index to compute thresholds
    !  - A split after xs(i) creates nodes of size i and (n_obs - i).
    !  - For a valid split, we need i >= min_obs and (n_obs - i) >= min_obs.
    !   This means we only need to consider unique values between xs(min_obs)
    !   and xs(n_train - min_obs + 1).
    n_start = min_obs                ! first X to use
    n_consider = count(x_mask)       ! number of observations to consider
    n_end = n_consider - n_start + 1 ! last X
    thresholds%nt = 0

    ! Early exit if a split is impossible from the start.
    if (n_end <= n_start) return

    ! Copy the relevant X values and sort them
    xs = pack(x(:, f), x_mask)
    xs = sort(n_consider, xs)

    ! loop over all feature values.
    !  - Skip over adjacent feature values that are too close
    !  - If a new unique value is found store it.
    n_unique = 1
    xs(1) = xs(n_start)
    do i = n_start + 1, n_end
      if (abs(xs(i) - xs(n_unique)) < feature_threshold) cycle
      n_unique = n_unique + 1
      xs(n_unique) = xs(i)
    end do

    ! If any threshould was found compute the thresholds
    !    0.5 * (Xs(i) + Xs(i + 1))
    ! and save thr (allocatable component in the type info_thr structure)
    thresholds%nt = n_unique - 1
    if (thresholds%nt > 0) thresholds%thr = 0.5_dp * (xs(1:n_unique - 1) + xs(2:n_unique))
  end subroutine find_thresholds

  pure subroutine select_imp_features(fa_id, feature, n_feat, bounds, imp_feat, n_imp)
    !-----------------------------------------------------------------------------
    ! Find the important features to avoid extra calculations
    !  - current feature is always important
    !  - ignore non missing feature if X_f is the range is (-Inf, Inf)
    !
    ! ARGUMENTS
    !   fa_id        [in]  : Integer, id of the father node
    !   feature      [in]  : Integer, index of the current feature
    !   n_feat       [in]  : Integer, number of features
    !   bounds       [in]  : Real(dp), entire bounds matrix
    !   imp_feat     [out] : Integer(n_feat), fixed-size array of important features
    !   n_imp        [out] : Integer, number of important features found
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: fa_id, feature, n_feat
    real(dp), intent(in) :: bounds(:, :)
    integer, intent(out) :: imp_feat(:)
    integer, intent(out) :: n_imp
    integer :: f, offset

    ! Earlie exit if there is only one feature
    if (n_feat == 1) then
      n_imp = 1
      imp_feat(1) = 1
      return
    end if

    ! Find the index of important features
    ! A feature is important if it's the current split feature, or if
    ! it has been used in a previous split (i.e., has finite bounds).
    offset = (fa_id - 1) * n_feat
    n_imp = 0
    do f = 1, n_feat
      if (f == feature .or. bounds(offset + f, 1) > neg_inf .or. bounds(offset + f, 2) < pos_inf) then
        n_imp = n_imp + 1
        imp_feat(n_imp) = f
      end if
    end do
  end subroutine select_imp_features

  pure subroutine update_imp_features(feature, n_imp, imp_feat)
    !-------------------------------------------------------------------------------
    ! Update the important features vector adding feature f.
    !
    ! ARGUMENTS
    !   feature      [in]  : Integer, index of the current feature
    !   n_imp        [out] : Integer, number of important features previously found
    !   imp_feat   [inout] : Integer, pre-allocated array of important features
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: feature
    integer, intent(inout) :: imp_feat(:)
    integer, intent(inout) :: n_imp
    integer :: f

    ! Earlie exit if there was no important feature
    if (n_imp == 0) then
      imp_feat = 0
      n_imp = 1
      imp_feat(1) = feature
      return
    end if

    ! Find the index of the first element greater than the new feature.
    ! If no such element exists, f will be n_imp + 1.
    f = count(imp_feat(1:n_imp) < feature) + 1

    ! Shift elements to the right in-place and insert the new feature
    if (f <= n_imp) imp_feat(f + 1 : n_imp + 1) = imp_feat(f : n_imp)
    imp_feat(f) = feature
    n_imp = n_imp + 1
  end subroutine update_imp_features

  function pdist(argsd, n, v, mu, sigma) result(fn_val)
    !------------------------------------------------------------------------
    ! Dispatches to the correct CDF function based on dist_id.
		!
		! Uses the fact that
		!     V = Z * sigma + mu,   Z ~ P_par
		! so that
		!     P(V <= v) = P_par(Z <= (v - mu)/sigma)
    !
    ! ARGUMENTS:
    !   argsd [in] : argsDist, distribution related parameters.
    !   n     [in] : Integer, number of observations.
    !   v     [in] : Real(dp), argument to compute the probability P(V <= v).
    !   mu    [in] : Real(dp), distribution parameters.
    !  sigma  [in] : Real(dp), distribution parameters.
    !
    ! Added July/2025
		!------------------------------------------------------------------------
		! Last revision: February, 2026
		!  - corrected the scaling: old version has v * sigma + mu
		!     instead of (v - mu)/sigma
    !------------------------------------------------------------------------
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n
    real(dp), intent(in) :: v, mu(n), sigma
    real(dp), allocatable :: fn_val(:)

    allocate(fn_val(n))

    ! Imitate CART (deterministic step function) when sigma is zero or near zero
    if (sigma <= eps) then
      where (mu <= v)
        fn_val = 1.0_dp
      elsewhere
        fn_val = 0.0_dp
      end where
      return
    end if

    select case (argsd%dist_id)
    case (1) ! norm
      ! P_par = N(0,1) - No need to transform
      ! Computes the cumulative probability P(V <= v) for V ~ N(mu, sigma^2).
      call pnorm_v(n, v, mu, sigma, fn_val)
    case (2) ! lnorm
      ! Computes the cumulative probability P(Z <= z)
      ! where z = (v - mu) * inv_sigma
      call plnorm_v(n, v, mu, sigma, argsd%dist_par, fn_val)
    case (3) ! t
      ! Computes the cumulative probability P(Z <= z)
      ! where z = (v - mu) * inv_sigma
      call pt_v(n, v, mu, sigma, argsd%dist_par, fn_val)
    case (4) ! gamma
      ! Computes the cumulative probability P(Z <= z)
      ! where z = (v - mu) * inv_sigma
      call pgamma_v(n, v, mu, sigma, argsd%dist_par, fn_val)
    end select
  end function pdist
end module prtree_misc
