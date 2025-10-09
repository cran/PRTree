subroutine pr_tree_fort(n_obs, n_feat, n_train, idx_train, y, x, n_sigmas, sigmas, int_param, dble_param, &
                        n_tn, p, gammahat, yhat, mse, nodes_info, thresholds, bounds, sigma_best, xregion)
  !------------------------------------------------------------------------------------
  ! Probabilistic Regression Tree
  !
  ! ARGUMENTS
  !   n_obs      [in]  : Integer, number of observations.
  !   n_feat     [in]  : Integer, number of features.
  !   n_train    [in]  : Integer, number of training observations.
  !   idx_train  [in]  : Integer, indices of the training observations.
  !   y          [in]  : Real(dp), response vector.
  !   X          [in]  : Real(dp), feature matrix.
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
  !   dble_param [in]   : Real(dp), parameters to build the tree
  !                     - min_prop: minimum proportion for splitting.
  !                     - min_prob: minimum probability for splitting.
  !                     - cp: complexity parameter (stopping criterion).
  !                     - my_par: parameters for the distribution.
  !   n_tn       [out] : Integer, number of final regions in the tree.
  !   P          [out] : Real(dp), probability matrix for the final tree (n_obs x max_tn).
  !   gammahat   [out] : Real(dp), vector of coefficients for the final tree (max_tn).
  !   yhat       [out] : Real(dp), predicted values for Y (n_obs).
  !   mse        [out] : Real(dp), mean square error for the final tree.
  !   nodes_info [out] : Integer, matrix with information about the nodes (2*max_tn-1, 5).
  !   thresholds [out] : Real(dp), vector of thresholds defining the regions (2*max_tn-1).
  !   bounds     [out] : Real(dp), lower and upper for regions (n_feat*(2*max_tn-1), 2).
  !   sigma_best [out] : Real(dp), selected sigma value in the grid (n_feat).
  !   Xregion    [out] : Integer, final region assigned to each X (n_obs).
  !
  ! DETAILS
  !  - For each sigma in the grid, builds a tree and selects the one with the lowest MSE.
  !  - When there are missing values, a proxy is used to assign X to a region.
  !  - The kth reg will be considered for splitting if, and only if,
  !    at least min_prop rows satisfy P(i,k) > min_prob.
  !  - The algorithm stops if 1 - mse_new/mse_old < cp.
  !
  ! HISTORY
  !  July/2025
  !     - Added na, which is then passed to main_calc and used by other functions.
  !     - var_inf and var_sup were replaced by bounds.
  !     - Added n_cand.
  !------------------------------------------------------------------------------------
  use prtree_main
  implicit none
  ! Input
  integer :: int_param(9)
  integer :: n_obs, n_feat, n_train, n_sigmas
  integer :: idx_train(n_train)
  real(dp) :: y(n_obs), x(n_obs, n_feat), sigmas(n_sigmas, n_feat)
  real(dp) :: dble_param(4)
  ! Output
  integer :: n_tn
  real(dp) :: p(n_obs, int_param(3)), mse(3)
  real(dp) :: gammahat(int_param(3)), yhat(n_obs)
  integer :: nodes_info(2 * int_param(3) - 1, 5)
  real(dp) :: thresholds(2 * int_param(3) - 1)
  real(dp) :: bounds(n_feat * (2 * int_param(3) - 1), 2)
  real(dp) :: sigma_best(n_feat)
  integer :: xregion(n_obs)
  ! Auxiliar variables
  real(dp) :: best_mse
  integer :: i, idx_test(max(1, n_obs - n_train)), n_test
  type(tree_model) :: tree
  type(tree_net) :: net_best, net_base
  type(tree_tdata) :: test_best

  ! print info level
  printinfo = int_param(9)

  ! find indexes for the test set
  idx_test = get_idx_test(n_obs, n_train, idx_train)
  n_test = n_obs - n_train

  ! Initializing the tree model
  call set_tree_model(tree, n_train, n_feat, y(idx_train), x(idx_train, :), int_param, dble_param)

  if (printinfo >= log_base) then
    ! Debug: Print initialization parameters
    call labelpr("===================================== ", -1)
    call labelpr("===  Initializing PRTree model == ", -1)
    call labelpr("===================================== ", -1)
    call intpr1("    Observations (n_obs):", -1, tree%train%n_obs)
    call intpr1("    Features (n_feat):", -1, tree%train%n_feat)
    call intpr1("    Max terminal nodes (max_tn):", -1, tree%ctrl%max_tn)
    call intpr1("    Min observations (min_obs):", -1, tree%ctrl%min_obs)
    call dblepr1("    Min proportion (min_prop):", -1, tree%ctrl%min_prop)
    call dblepr1("    Min probability (min_prob):", -1, tree%ctrl%min_prob)
    call dblepr1("    Complexity parameter (cp):", -1, tree%ctrl%cp)
    call intpr1("     Filling type (fill):", -1, tree%train%fill)
    call intpr1("     Proxy criterion (crit: 1 = mean, 2 = variance, 3 = both):", -1, tree%ctrl%crit)
    call labelpr(" ", -1)

    call labelpr("==============================================", -1)
    call labelpr("===== Setting the distribution parameters ====", -1)
    call labelpr("==============================================", -1)
    select case (int_param(8))
    case (1)
      call labelpr("Using the standard Gaussin distribution", -1)
    case (2)
      call dblepr1("Using the Log-normal distribution with meanlog = 0 and sdlog = ", -1, tree%dist%dist_par)
    case (3)
      call dblepr1("Using the Student's t distribution with df = ", -1, tree%dist%dist_par)
    case (4)
      call dblepr1("Using the Gamma distribution with scale = 1 and shape = ", -1, tree%dist%dist_par)
    end select
    call labelpr(" ", -1)
  end if

  if (n_test > 0) then
    ! Initializes the validation data set
    tree%test%n_obs = n_test
    tree%test%x = x(idx_test, :)
    tree%test%y = y(idx_test)
    tree%test%yhat = vector(n_test, source=0.0_dp)
  end if

  ! loop to select the search sigma
  best_mse = pos_inf

  ! Save the initial status to avoid recomputing the thresholds for each sigma
  call start_root(tree)

  ! Early skip if no thresholds are available
  if (all(tree%net%nodes(1)%state%thresholds%nt == 0)) then
    ! Debug: Print node splitting info
    call labelpr("Root node cannot be splitted: No thresholds found.", -1)
    call labelpr(" ", -1)

    ! The tree is just the root node. Populate outputs accordingly before returning.
    sigma_best = 0.0_dp

    ! Populate output arrays for the training set based on the single-node tree
    call return_tree(tree%net, n_obs, tree%train%n_feat, n_train, idx_train, tree%ctrl%max_tn, n_tn, &
                     nodes_info, thresholds, bounds, p, gammahat, yhat, mse(1), xregion)

    ! Handle test set if it exists
    if (n_test > 0) then
      ! For a single-node tree, prediction is the mean of the training data
      yhat(idx_test) = gammahat(1)
      mse(2) = sum((y(idx_test) - gammahat(1))**2)
      mse(3) = (tree%net%mse + mse(2)) / n_obs
      mse(2) = mse(2) / n_test
      p(idx_test, 1) = 1.0_dp
      xregion(idx_test) = 0
    else
      mse(2) = 0.0_dp
      mse(3) = mse(1)
    end if
    return
  end if

  if (printinfo >= log_base) then
    ! Debug: Sigma grid search
    call labelpr("====================================", -1)
    call labelpr("==== Starting sigma grid search ====", -1)
    call labelpr("====================================", -1)
    call labelpr(" ", -1)
  end if

  net_base = tree%net
  do i = 1, n_sigmas
    ! fix sigma and build the candidate tree
    tree%train%sigma = sigmas(i, :)

    if (printinfo >= log_verbose) then
      ! Debug: Sigma grid search
      call dblepr("Trying sigma =", -1, tree%train%sigma, n_feat)
      call labelpr(" ", -1)
    end if

    call build_tree(tree)
    if (n_test > 0) then
      ! compute the predicitons and the mse for the validation/testing sample
      call predict_tree(tree)
    else
      ! use the mse from the training set
      tree%test%mse = tree%net%mse
    end if

    if (printinfo >= log_verbose .and. i > 1) then
      ! Debug: Sigma grid search
      call dblepr1("Current MSE:", -1, best_mse)
      call dblepr1("Competing MSE:", -1, tree%test%mse)
      call labelpr(" ", -1)
    end if

    ! For each sigma fixed, the search trees is obtained.
    ! Updates the current tree when a better one is found.
    if (i == 1 .or. tree%test%mse < best_mse) then
      if (printinfo >= log_base) then
        ! Debug: Sigma grid search
        if (i == 1) then
          call labelpr("First tree status.", -1)
        else
          call labelpr("New best sigma found.", -1)
        end if
        call dblepr("    sigma ", -1, tree%train%sigma, n_feat)
        call dblepr1("    Best sigma (proxy) MSE (test):", -1, tree%test%mse)
        call labelpr(" ", -1)
      end if
      best_mse = tree%test%mse
      sigma_best = tree%train%sigma
      net_best = tree%net
      if (n_test > 0) test_best = tree%test
    end if
    if (n_sigmas > 1) tree%net = net_base
  end do

  ! updating values and returning to R (training and test sets)
  call return_tree(net_best, n_obs, tree%train%n_feat, n_train, idx_train, tree%ctrl%max_tn, n_tn, &
                   nodes_info, thresholds, bounds, p, gammahat, yhat, mse(1), xregion)

  if (n_test > 0) then
    call return_tree(test_best, n_obs, n_test, idx_test, n_tn, p(1:n_obs, 1:n_tn), yhat, mse(2), xregion)

    mse(3) = (net_best%mse + test_best%mse) / n_obs
  else
    mse(2) = 0.0_dp
    mse(3) = mse(1)
  end if

  if (printinfo >= log_base) then
    ! Debug: Final tree summary
    call labelpr("=== PRTree construction complete.==  ", -1)
    call intpr1("    Terminal nodes (n_tn):", -1, n_tn)
    call dblepr1("    Final MSE (training set):", -1, mse(1))
    if (n_test > 0) then
      call dblepr1("    Final MSE (test set):", -1, mse(2))
      call dblepr1("    Final MSE (global):", -1, mse(3))
    end if
    call intpr("    Terminal node IDs:", -1, net_best%tn_id, n_tn)
  end if
end subroutine pr_tree_fort

subroutine predict_pr_tree_fort(dist, pardist, fill, n_obs, n_feat, x_test, bounds, n_tn, tn, nodes_info, &
                                p, gammahat, sigma, yhat_test)
  !--------------------------------------------------------------------------------------
  ! Predicts the response for new data using the probabilistic regression tree.
  !
  ! ARGUMENTS
  !   dist       [in]  : Integer, the code for the distribution
  !   pardist    [in]  : Real(dp), parameters for the distribution (if any)
  !   fill       [in]  : Integer, method used to fill P(i,:)
  !                      when there is at least one missing feature in X(i,:).
  !   n_obs      [in]  : Integer, number of observations in the test
  !   n_feat     [in]  : Integer, number of features.
  !   X_test     [in]  : Real(dp), feature matrix for the test data (n_obs x n_feat).
  !   bounds     [in]  : Real(dp), bounds for each feature in the nodes
  !                      (n_feat*(2*n_tn-1), 2).
  !   n_tn       [in]  : Integer, number of terminal nodes in the tree.
  !   tn         [in]  : Integer, terminal node indices (n_tn).
  !   nodes_info [in]  : Integer, information on nodes (2*n_tn-1).
  !   P          [out] : Real(dp), probability matrix for the test data (n_obs x n_tn).
  !   gammahat   [in]  : Real(dp), vector of coefficients for the final tree (n_tn).
  !   sigma      [in]  : Real(dp), vector of standard deviations for the features (n_feat).
  !   yhat_test  [out] : Real(dp), predicted values for Y in the test data (n_obs).
  !
  !   DETAILS
  !     - For each observation, computes the probabilities of belonging to each terminal node.
  !     - Uses the bounds to compute the probabilities.
  !     - Fills the matrix P with the probabilities.
  !--------------------------------------------------------------------------------------
  use prtree_main
  implicit none
  ! Input
  integer :: dist, n_obs, n_feat, n_tn, fill
  integer :: tn(n_tn), nodes_info(2 * n_tn - 1, 4)
  real(dp) :: x_test(n_obs, n_feat), gammahat(n_tn)
  real(dp) :: sigma(n_feat), pardist
  real(dp) :: bounds(n_feat * (2 * n_tn - 1), 2)
  ! Output
  real(dp) :: yhat_test(n_obs)
  real(dp) :: p(n_obs, n_tn)
  ! Local
  type(argsdist) :: argsd
  !---------------------------------------------
  ! filling the matrix P with the values
  !    Psi(x, Rj, sigma), 1 <= j <= n_tn
  ! where Rj is the j-th terminal node (reg)
  !---------------------------------------------
  argsd%dist_id = dist
  argsd%dist_par = pardist
  p = prob_rgivenx(argsd, n_obs, n_feat, x_test, sigma, fill, n_tn, tn, nodes_info, bounds)
  yhat_test = matmul(p, gammahat)
end subroutine predict_pr_tree_fort
