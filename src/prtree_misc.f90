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
    real(dp) :: xout(n)
    integer :: i, idx(n)

    ! make a copy to avoid changing the original data
    xout = x
    if (n < 2) return

    ! Use the more efficient Heapsort (revsortr) and reverse the result.
    idx = [(i, i=1, n)]
    call revsortr(xout, idx, n)  ! Sorts descending

    ! Reverse in-place to get ascending order
    xout = xout(n:1:-1)
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
    integer :: idx_n(n)

    if (n < 2) return

    ! 1. Create an index array [1, 2, ..., n]
    idx_n = [(i, i=1, n)]

    ! 2. Sort x in descending order and permut idx to match.
    call revsortr(x, idx_n, n)

    ! 3. Reorder the original y and idx values according to the permuted index.
    y = y(idx_n)
    idx = idx(idx_n)

    ! 4. Reverse the arrays in-place to get the final ascending order.
    x = x(n:1:-1)
    y = y(n:1:-1)
    idx = idx(n:1:-1)
  end subroutine sort_xy

  function lsquare(m, n, a, b) result(x)
    !-------------------------------------------------------------------------------------
    ! Use LAPACK's DGELSD to solve the least squares problem
    !   minimum || b - A * x||
    !
    ! ARGUMENTS
    !   m [in] : Integer, number of rows
    !   n [in] : Integer, number of columns
    !   A [in] : Real(dp), m x n matrix
    !   b [in] : Real(dp), size m vector
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
    integer, intent(in) :: m, n
    real(dp), intent(in) :: a(m, n), b(m)
    real(dp) :: x(n)
    real(dp) :: a_copy(m, n), b_copy(max(m, n), 1)
    real(dp) :: rcond, s(min(m, n))
    integer :: lda, ldb, info, rank, lwork
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    rcond = -1.0_dp  ! rcond < 0 means: use machine precision
    lda = m
    ldb = max(m, n)

    ! Copy input data (DGELSD overwrites its input arrays)
    a_copy = a
    b_copy(1:m, 1) = b

    ! First call: workspace query to find the optimum working size
    lwork = -1
    info = 0
    allocate(work(1)) ! to avoid error when compiling for mac
    work = 0.0_dp
    allocate(iwork(1)) ! to avoid error when compiling for mac
    iwork = 1
    call dgelsd(m, n, 1, a_copy, lda, b_copy, ldb, s, rcond, rank, work, lwork, iwork, info)

    if (info /= 0) then
      x = 0.0_dp
      return
    end if

    ! reallocate work and solve the least square problem
    lwork = max(1, int(work(1)))
    work = vector(lwork, 0.0_dp)
    iwork = vector(iwork(1), 0)
    call dgelsd(m, n, 1, a_copy, lda, b_copy, ldb, s, rcond, rank, work, lwork, iwork, info)

    ! The solution is returned in the first k rows of B
    if (info == 0) then
      x = b_copy(1:n, 1)
    else
      x = 0.0_dp
      if (printinfo >= log_base) then
        call labelpr("Warning: DGELSD failed in `lsquare`", -1)
        call intpr1("Info code:", -1, info)
      end if
    end if
  end function lsquare

  pure function get_idx_test(n_obs, n_train, idx_train) result(idx_test)
    !--------------------------------------------------------------------------
    !  Returns the indexes of the testing set
    !
    ! ARGUMENTS
    !   n_obs     [in] : Integer, total number of observations.
    !   n_train   [in] : Integer, number of observations in the training set.
    !   idx_train [in] : Integer, array of indexes for the training set.
    !
    ! RETURNS
    !   idx_test : Integer, array of indexes for the testing set.
    !
    ! Added July/2025
    !--------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n_obs, n_train
    integer, intent(in) :: idx_train(n_train)
    integer :: idx_test(max(1, n_obs - n_train))
    logical :: is_train_mask(n_obs)
    integer :: i

    ! if n_train = n_obs, no test set is needed
    if (n_train == n_obs) then
      idx_test = 0
      return
    end if

    ! mask for indexes in the training set
    is_train_mask = .false.
    is_train_mask(idx_train) = .true.
    idx_test = pack([(i, i=1, n_obs)], mask=.not. is_train_mask)
  end function get_idx_test

  subroutine find_thresholds(min_obs, n_obs, x, thresholds)
    !------------------------------------------------------------------------------
    ! Finds candidate thresholds for splitting a node. A threshold is the midpoint
    ! between two unique, adjacent feature values. A split is only valid if it
    ! results in child nodes each having at least min_obs observations.
    !
    ! ARGUMENTS
    !   min_obs    [in]  : Integer, minimum observations per node.
    !   n_obs      [in]  : Integer, number of observations.
    !   x          [in]  : Real(dp), feature values.
    !   thresholds [out] : Type(info_thr), thresholds found for the feature.
    !
    ! Added July/2025
    !------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n_obs, min_obs
    real(dp), intent(in) :: x(n_obs)
    type(info_thr), intent(out) :: thresholds
    real(dp) :: xs(n_obs)
    integer :: i, n_start, n_end, n_unique

    ! Find the first and last index to compute thresholds
    !  - A split after xs(i) creates nodes of size i and (n_obs - i).
    !  - For a valid split, we need i >= min_obs and (n_obs - i) >= min_obs.
    !   This means we only need to consider unique values between xs(min_obs)
    !   and xs(n_obs - min_obs + 1).
    n_start = min_obs           ! first X
    n_end = n_obs - n_start + 1 ! last X
    thresholds%nt = 0

    ! Early exit if a split is impossible from the start.
    if (n_end <= n_start) return

    ! Sort the X values
    xs = sort(n_obs, x)

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

  pure subroutine select_imp_features(feature, n_feat, bounds, imp_feat, n_imp)
    !-----------------------------------------------------------------------------
    ! Find the important features to avoid extra calculations
    !  - current feature is always important
    !  - ignore non missing feature if X_f is the range is (-Inf, Inf)
    !
    ! ARGUMENTS
    !   feature      [in]  : Integer, index of the current feature
    !   n_feat       [in]  : Integer, number of features
    !   bounds       [in]  : Real(dp), bounds for each feature (n_feat x 2)
    !   imp_feat     [out] : Integer, array of important features
    !   n_imp        [out] : Integer, number of important features found
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: feature, n_feat
    real(dp), intent(in) :: bounds(n_feat, 2)
    integer, allocatable, intent(out) :: imp_feat(:)
    integer, intent(out) :: n_imp
    integer :: f
    integer :: temp(n_feat)

    ! Earlie exit if there is only one feature
    if (n_feat == 1) then
      n_imp = 1
      imp_feat = [1]
      return
    end if

    ! Find the index of important features
    ! A feature is important if it's the current split feature, or if
    ! it has been used in a previous split (i.e., has finite bounds).
    n_imp = 0
    do f = 1, n_feat
      if (f == feature .or. bounds(f, 1) > neg_inf .or. bounds(f, 2) < pos_inf) then
        n_imp = n_imp + 1
        temp(n_imp) = f
      end if
    end do

    ! Copy the collected indices to the final allocatable array.
    if (n_imp > 0) imp_feat = temp(1:n_imp)
  end subroutine select_imp_features

  pure subroutine update_imp_features(feature, n_imp, imp_feat)
    !-------------------------------------------------------------------------------
    ! Update the important features vector adding feature f.
    !
    ! ARGUMENTS
    !   feature      [in]  : Integer, index of the current feature
    !   n_imp        [out] : Integer, number of important features previously found
    !   imp_feat     [out] : Integer, array of important features
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: feature
    integer, allocatable, intent(inout) :: imp_feat(:)
    integer, intent(inout) :: n_imp
    integer :: f
    integer :: temp(n_imp + 1)

    ! Earlie exit if there was no important feature
    if (n_imp == 0) then
      n_imp = 1
      imp_feat = [feature]
      return
    end if

    ! Find the index of the first element greater than the new feature.
    ! If no such element exists, f will be n_imp + 1.
    f = count(imp_feat < feature) + 1

    ! Construct the new array by inserting the feature at the correct position.
    temp(1 : f - 1) = imp_feat(1 : f - 1)
    temp(f) = feature
    temp(f + 1 : n_imp + 1) = imp_feat(f : n_imp)

    ! Update the original array and its size.
    n_imp = n_imp + 1
    imp_feat = temp
  end subroutine update_imp_features

  pure function pdist(argsd, n, x, mu, sigma) result(fn_val)
    !------------------------------------------------------------------------
    ! Dispatches to the correct CDF function based on dist_id.
    !
    ! ARGUMENTS:
    !   argsd [in] : argsDist, distribution related parameters.
    !   n     [in] : Integer, number of observations.
    !   x     [in] : Real(dp), argument to compute the probability P(X <= x).
    !   mu    [in] : Real(dp), distribution parameters.
    !  sigma  [in] : Real(dp), distribution parameters.
    !
    ! Added July/2025
    !------------------------------------------------------------------------
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n
    real(dp), intent(in) :: x, mu(n), sigma
    real(dp) :: fn_val(n), xsig
    integer :: i

    ! Precompute transformation if needed
    if (argsd%dist_id /= 1) xsig = x * sigma

    select case (argsd%dist_id)
    case (1) ! norm
      ! Computes the cumulative probability P(X <= x) for X ~ N(mu, sigma^2).
      do i = 1, n
        fn_val(i) = pnorm_pure(x, mu(i), sigma)
      end do
    case (2) ! lnorm
      ! Computes the cumulative probability P(X <= x)
      ! where X = (Z - mu)/sigma and Z ~ lognormal(meanlog = 0, sdlog = sdlog).
      do i = 1, n
        fn_val(i) = plnorm_pure(xsig + mu(i), 0.0_dp, argsd%dist_par)
      end do
    case (3) ! t
      ! Computes the cumulative probability P(X <= x)
      ! where X = (Z - mu)/sigma and X ~ t(df = df).
      do i = 1, n
        fn_val(i) = pt_pure(xsig + mu(i), argsd%dist_par)
      end do
    case (4) ! gamma
      ! Computes the cumulative probability P(X <= x)
      ! where X = (Z - mu)/sigma and Z ~ Gamma(shape, scale = 1).
      do i = 1, n
        fn_val(i) = pgamma_pure(xsig + mu(i), argsd%dist_par, 1.0_dp)
      end do
    end select
  end function pdist
end module prtree_misc
