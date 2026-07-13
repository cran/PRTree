subroutine pr_tree_fort(n_obs, n_feat, n_train, y_train, y_test, x_train, x_test, &
  n_sigmas, sigmas, int_param, dble_param, n_tn, p_train, p_test, gammahat, yhat_train, &
  yhat_test, mse, nodes_info, thresholds, sigma_best, xregion)
  !------------------------------------------------------------------------------------
  ! Probabilistic Regression Tree
  !
  ! ARGUMENTS
  !   n_obs      [in]  : Integer, number of observations.
  !   n_feat     [in]  : Integer, number of features.
  !   n_train    [in]  : Integer, number of training observations.
  !   y_train    [in]  : Real(dp), response vector (train set).
  !   y_test     [in]  : Real(dp), response vector (test set).
  !   x_train    [in]  : Real(dp), feature matrix (train set).
  !   x_test     [in]  : Real(dp), feature matrix (test set).
  !   n_sigmas   [in]  : Integer, number of sigma values in the grid.
  !   sigmas     [in]  : Real(dp), matrix (n_sigmas x n_feat) of sigma values.
  !   int_param  [in]  : integer(:), parameters to build the tree
  !                     - fill: method used to fill P(i,:) when there are missing features
  !                     - crit: criterion used to assign a missing value to a region.
  !                     - max_tn: maximum number of terminal nodes allowed.
  !                     - max_d: maximum tree depth.
  !                     - min_obs: minimum number of observations in a final node.
  !                     - n_cand: number of candidate thresholds for split search.
  !                     - by_node: if 1 then the best candidates are found by node.
  !                     - my_dist: the code for the distribution.
  !                     - iprint: controls the print detail level.
  !                     - output_mode: 0 = CV/light output, 1 = full model
  !   dble_param [in]  : Real(dp), parameters to build the tree
  !                     - min_prop: minimum proportion for splitting.
  !                     - min_prob: minimum probability for splitting.
  !                     - cp: complexity parameter (stopping criterion).
  !                     - my_par: parameters for the distribution.
  !   n_tn       [out] : Integer, number of final regions in the tree.
  !   P_train    [out] : Real(dp), probability matrix for the final tree (n_train x max_tn).
  !   P_test     [out] : Real(dp), probability matrix for the final tree (n_test x max_tn).
  !   gammahat   [out] : Real(dp), vector of coefficients for the final tree (max_tn).
  !   yhat_train [out] : Real(dp), predicted values for Y (n_train).
  !   yhat_test  [out] : Real(dp), predicted values for Y (n_test).
  !   mse        [out] : Real(dp), mean square error for the final tree.
  !   nodes_info [out] : Integer, matrix with information about the nodes (2*max_tn-1, 5).
  !   thresholds [out] : Real(dp), vector of thresholds defining the regions (2*max_tn-1).
  !   sigma_best [out] : Real(dp), selected sigma value in the grid (n_feat).
  !   Xregion    [out] : Integer, final region assigned to each X (n_train).
  !
  ! DETAILS
  !  - For each sigma in the grid, builds a tree and selects the one with the lowest MSE.
  !  - When there are missing values, a proxy is used to assign X to a region.
  !  - The kth reg will be considered for splitting if, and only if,
  !    at least min_prop rows satisfy P(i,k) > min_prob.
  !  - The algorithm stops if 1 - mse_new/mse_old < cp.
  !
  ! MEMORY DESIGN
  !  - No base copy of tree is stored between sigma values.
  !  - The root state is rebuilt with start_root(tree) before each new sigma.
  !  - In full-output mode, the best tree is serialized immediately to the
  !    output buffers already allocated by R.
  !
  ! HISTORY
  !  July/2025
  !     - Added na, which is then passed to main_calc and used by other functions.
  !     - var_inf and var_sup were replaced by bounds.
  !     - Added n_cand.
  !  March/2026
  !     - Refactored to avoid keeping a full copy of tree in memory
  !       during the sigma search.
  !  May/2026
  !    - Created pointers to reduce the number of copies and reduce memory usage.
  !    - Now train and testing arguments are passed separately instead of as a single
  !      matrix/vector (needed for pointers).
  !    - Removed the bounds argument to use less memory. The bounds are now computed
  !      on the fly when needed.
  !------------------------------------------------------------------------------------
  use prtree_main
  implicit none

  ! Input
  integer :: int_param(10)
  integer :: n_obs, n_feat, n_train, n_sigmas
  real(dp), target :: y_train(n_train)
  real(dp), target :: y_test(max(1, n_obs - n_train))
  real(dp), target :: x_train(n_train, n_feat)
  real(dp), target :: x_test(max(1, n_obs - n_train), n_feat)
  real(dp) :: sigmas(n_sigmas, n_feat)
  real(dp) :: dble_param(4)

  ! Output
  integer :: n_tn
  real(dp) :: mse(3)
  real(dp) :: sigma_best(n_feat)
  integer,  target :: xregion(max(1, n_train * int_param(10)))
  real(dp), target :: gammahat(max(1, int_param(3) * int_param(10)))
  real(dp), target :: yhat_train(max(1, n_train * int_param(10)))
  real(dp), target :: yhat_test(max(1, (n_obs - n_train) * int_param(10)))
  real(dp), target :: thresholds(max(1, (2 * int_param(3) - 1) * int_param(10)))
  real(dp), target :: p_train(max(1, n_train * int_param(10)), &
                              max(1, int_param(3) * int_param(10)))
  real(dp), target :: p_test(max(1, (n_obs - n_train) * int_param(10)), &
                             max(1, int_param(3) * int_param(10)))
  integer,  target :: nodes_info(max(1, (2 * int_param(3) - 1) * int_param(10)), &
                                 max(1, 5 * int_param(10)))

  ! Auxiliar variables
  real(dp) :: best_test_mse, best_train_mse
  integer :: i
  type(tree_model), target :: tree
  logical :: keep_full

  ! set tree parameters that do not depend on sigma
  call set_tree_model(tree, int_param, dble_param, keep_full, n_obs, n_train, n_feat, &
    y_train, x_train, y_test, x_test)
  call assign_pointers(tree, keep_full, n_sigmas, yhat_train, yhat_test, xregion, &
     p_train, p_test, gammahat, nodes_info, thresholds)

  ! Build the root structure for the first sigma
  call start_root(tree)

  ! Initialize best metrics
  best_test_mse = pos_inf
  ! Early skip if no thresholds are available
  if (all(tree%nodes(1)%thresholds%nt == 0)) then
    call return_root_only_tree(tree, y_test, keep_full, n_tn, yhat_train, yhat_test, &
       p_test, mse, nodes_info, thresholds, sigma_best, xregion)
    return
  end if

   if (tree%printinfo >= log_base) then
    ! Debug: Sigma grid search
    call labelpr("====================================", -1)
    call labelpr("==== Starting sigma grid search ====", -1)
    call labelpr("====================================", -1)
    call labelpr(" ", -1)
  end if

  do i = 1, n_sigmas
    ! Rebuild the root state before each new sigma after the first one.
    ! This avoids keeping a full base copy of tree in memory.
    if (i > 1) call start_root(tree)

    ! Fix sigma and build the candidate tree
    tree%sigma = sigmas(i, :)

    if (tree%printinfo >= log_verbose) then
      ! Debug: Sigma grid search
      call dblepr("Trying sigma =", -1, tree%sigma, n_feat)
      call labelpr(" ", -1)
    end if

    ! build the tree for the current sigma
    ! predict: if there is no test set, just copy the mse
    call build_tree(tree)
    call predict_tree(tree)

    if (tree%printinfo >= log_verbose .and. i > 1) then
      ! Debug: Sigma grid search
      call dblepr1("Current best MSE:", -1, best_test_mse)
      call dblepr1("Competing MSE:", -1, tree%mse_test)
      call labelpr(" ", -1)
    end if

    ! Update the best sigma whenever the current candidate improves the selection criterion.
    if (i == 1 .or. tree%mse_test < best_test_mse) then
      if (tree%printinfo >= log_base) then
        ! Debug: Sigma grid search
        if (i == 1) then
          call labelpr("First tree status.", -1)
        else
          call labelpr("New best sigma found.", -1)
        end if
        call dblepr("    sigma ", -1, tree%sigma, n_feat)
        call dblepr1("    Best sigma (proxy) MSE (test):", -1, tree%mse_test)
        call labelpr(" ", -1)
      end if

      ! Always store the scalar information required by CV/light mode.
      best_train_mse = tree%mse_train
      best_test_mse = tree%mse_test
      sigma_best = tree%sigma
      n_tn = tree%n_tn

      if (keep_full) then
        ! IMPORTANT:
        ! Instead of storing a deep copy of the tree, immediately serialize
        ! the current best tree to the output buffers that were already allocated by R.
        ! This keeps memory usage lower during the sigma search.
        call return_tree(tree, n_tn, nodes_info, thresholds, &
                         yhat_train, yhat_test, mse, xregion, i == n_sigmas)
      end if

    else if (i == n_sigmas .and. keep_full) then
      ! if n_sigmas = 1 the code will enter the previous if block and return the tree,
      ! so we only need to check keep_full here and, if we are not in light mode, we
      ! need to return the p_train, p_test and gammahat for the best tree found.
      p_train = 0.0_dp
      p_train(:, 1:n_tn) = tree%p_temp(:, 1:n_tn)
      if(tree%n_test > 0) then
        p_test = 0.0_dp
        p_test(:, 1:n_tn) = tree%p_test_temp(:, 1:n_tn)
      end if
      gammahat = 0.0_dp
      gammahat(1:n_tn) = tree%gamma_temp(1:n_tn)
    end if
  end do

  ! In light mode, no tree structure was serialized. Only return the metrics of the best sigma.
  if (.not. keep_full) then
    if (tree%n_test > 0) then
      mse(1) = best_train_mse / n_train
      mse(2) = best_test_mse / tree%n_test
      mse(3) = (best_train_mse + best_test_mse) / n_obs
    else
      mse(1) = best_test_mse / n_train
      mse(2) = 0.0_dp
      mse(3) = mse(1)
    end if
  end if

  if (tree%printinfo >= log_base) then
    ! Debug: Final tree summary
    call labelpr("=== PRTree construction complete.==  ", -1)
    call intpr1("    Terminal nodes (n_tn):", -1, n_tn)
    call dblepr1("    Final MSE (training set):", -1, mse(1))
    if (tree%n_test > 0) then
      call dblepr1("    Final MSE (test set):", -1, mse(2))
      call dblepr1("    Final MSE (global):", -1, mse(3))
    end if
  end if
end subroutine pr_tree_fort

subroutine predict_pr_tree_fort(dist, pardist, fill, n_test, n_feat, x_test, thresholds, &
                                n_tn, tn, nodes_info, p, gammahat, sigma, yhat_test)
  !--------------------------------------------------------------------------------------
  ! Predicts the response for new data using the probabilistic regression tree.
  !
  ! ARGUMENTS
  !   dist       [in]  : Integer, the code for the distribution
  !   pardist    [in]  : Real(dp), parameters for the distribution (if any)
  !   fill       [in]  : Integer, method used to fill P(i,:)
  !                      when there is at least one missing feature in X(i,:).
  !   n_test     [in]  : Integer, number of observations in the test
  !   n_feat     [in]  : Integer, number of features.
  !   X_test     [in]  : Real(dp), feature matrix for the test data (n_test x n_feat).
  !   thresholds [in]  : Real(dp), split thresholds for the tree nodes.
  !   n_tn       [in]  : Integer, number of terminal nodes in the tree.
  !   tn         [in]  : Integer, terminal node indices (n_tn).
  !   nodes_info [in]  : Integer, information on nodes (2*n_tn-1).
  !   P          [out] : Real(dp), probability matrix for the test data (n_test x n_tn).
  !   gammahat   [in]  : Real(dp), vector of coefficients for the final tree (n_tn).
  !   sigma      [in]  : Real(dp), vector of standard deviations for the features (n_feat).
  !   yhat_test  [out] : Real(dp), predicted values for Y in the test data (n_test).
  !
  !   DETAILS
  !     - For each observation, computes the probabilities of belonging to each terminal node.
  !     - Uses the bounds to compute the probabilities.
  !     - Fills the matrix P with the probabilities.
  !--------------------------------------------------------------------------------------
  use prtree_main
  implicit none
  ! Input
  integer :: dist, n_test, n_feat, n_tn, fill
  integer :: tn(n_tn), nodes_info(2 * n_tn - 1, 4)
  real(dp) :: x_test(n_test, n_feat), gammahat(n_tn)
  real(dp) :: sigma(n_feat), pardist
  real(dp) :: thresholds(2 * n_tn - 1)
  ! Output
  real(dp) :: yhat_test(n_test)
  real(dp) :: p(n_test, n_tn)
  ! Local
  type(argsdist) :: argsd
  !---------------------------------------------
  ! filling the matrix P with the values
  !    Psi(x, Rj, sigma), 1 <= j <= n_tn
  ! where Rj is the j-th terminal node (reg)
  !---------------------------------------------
  argsd%dist_id = dist
  argsd%dist_par = pardist
  p = prob_rgivenx(argsd, n_test, n_feat, x_test, sigma, fill, n_tn, tn, nodes_info, thresholds)
  yhat_test = matmul(p, gammahat)
end subroutine predict_pr_tree_fort
