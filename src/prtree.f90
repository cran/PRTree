!-------------------------------------------------------------------------------
! PRTree: Probabilistic Regression Tree
!-------------------------------------------------------------------------------
! Module prtree_types contains
!   - Derived types
!   - Interface for external functions
!   - Funtions to allocate/initialize vectors and matrices
!
! Module prtree_misc contains
!   - Utility functions
!
! Module prtree_main contains
!   - Building and reporting subroutines
!
! Description:
! Part 1: set module parameters, defines the types used in the PRTree algorithm,
! including the input data structure, control parameters, and tree structure.
! Also contains interfaces for external functions.
!
! Part 2: contains utility functions used in the PRTree algorithm, including
! sorting, sequence generation, vector/matrix initialization and least square
! estimation.
!
! Part 3: defines the main interface for the PRTree algorithm,
! including training, prediction, and model evaluation functions.
!
! Build upon the original PRTree implementation by Alisson S. Neimaier.
!
!-------------------------------------------------------------------------------
! This version: May/2026
! By Taiane S. Prass - PPGEst/UFRGS
!-------------------------------------------------------------------------------
module prtree_main
  use prtree_types
  use prtree_misc
  implicit none

contains
  pure subroutine get_feature_bounds(node_id, feat, nodes_info, thresholds, b_lower, b_upper)
    !------------------------------------------------------------------------------
    ! Reconstructs bounds only for a specific feature by tracing its ancestry.
    ! Eliminates the need for a global bounds matrix.
    !
    ! ARGUMENTS
    !   node_id    [in]  : Integer, id of the current node.
    !   feat       [in]  : Integer, id of the current feature.
    !   nodes_info [in]  : Integer, matrix with nodes information.
    !   thresholds [in]  : Real, vector of thresholds for all nodes.
    !   b_lower    [out] : Real, lower bound for the current feature.
    !   b_upper    [out] : Real, upper bound for the current feature.
    !
    ! NOTES
    !   - Left children are always EVEN, right children are always ODD.
    !
    ! Added May 2026
    !   - replaces the old bounds matrix with a more efficient approach that only
    !     retracks bounds for the current feature.
    !------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: node_id, feat
    integer, intent(in) :: nodes_info(:,:)
    real(dp), intent(in) :: thresholds(:)
    real(dp), intent(out) :: b_lower, b_upper
    integer :: curr, parent, split_feat, offset, ncol_ni

    b_lower = neg_inf
    b_upper = pos_inf
    ncol_ni = size(nodes_info, 2)
    offset = 5 - ncol_ni
    curr = node_id

    do while (curr > 1)
      parent = nodes_info(curr, idx_ni_father - offset)
      split_feat = nodes_info(parent, idx_ni_feature - offset)
      if (split_feat == feat) then
        ! Left children are always EVEN, right children are always ODD
        if (mod(curr, 2) == 0) then
          if (b_upper >= pos_inf) b_upper = thresholds(parent)
        else
          if (b_lower <= neg_inf) b_lower = thresholds(parent)
        end if
        if (b_lower > neg_inf .and. b_upper < pos_inf) exit
      end if
      curr = parent
    end do
  end subroutine get_feature_bounds

  pure subroutine get_node_imp_features(node_id, nodes_info, imp_feat, n_imp)
    !-------------------------------------------------------------------------------------
    ! Reconstructs important features list by tracing ancestry, replacing bounds array.
    ! Eliminates the need for a global bounds matrix.
    !
    ! ARGUMENTS
    !   node_id    [in]  : Integer, id of the current node.
    !   nodes_info [in]  : Integer, matrix with nodes information.
    !   imp_feat   [out] : Integer, vector of important features for the current node.
    !   n_imp      [out] : Integer, number of important features for the current node.
    !
    ! NOTES
    !   - Important features are collected by tracing the ancestry of the current node.
    !   - The list is unique and ordered by depth (closest ancestors first).
    !
    ! Added May 2026
    !  - replaces the old bounds matrix with a more efficient approach that only
    !  - retracks important features for the current node.
    !  - eliminates the need for a global important features matrix.
    !  - ensures that the important features list is unique and ordered by depth.
    !-------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: node_id
    integer, intent(in) :: nodes_info(:,:)
    integer, intent(out) :: imp_feat(:)
    integer, intent(out) :: n_imp
    integer :: curr, parent, split_feat, i, offset, ncol_ni
    logical :: already_added

    ncol_ni = size(nodes_info, 2)
    offset = 5 - ncol_ni
    n_imp = 0
    curr = node_id

    do while (curr > 1)
      parent = nodes_info(curr, idx_ni_father - offset)
      split_feat = nodes_info(parent, idx_ni_feature - offset)
      if (split_feat > 0) then
        already_added = .false.
        do i = 1, n_imp
          if (imp_feat(i) == split_feat) then
            already_added = .true.
            exit
          end if
        end do
        if (.not. already_added) then
          n_imp = n_imp + 1
          imp_feat(n_imp) = split_feat
        end if
      end if
      curr = parent
    end do
    if (n_imp == 0) then
      n_imp = 1
      imp_feat(1) = 1
    end if
  end subroutine get_node_imp_features

  !======================================================
  ! STEP 1: start the tree model
  ! Initialization performed:
  !======================================================
  subroutine print_starting_message(tree)
    !----------------------------------------------------------------------
    ! Prints the starting message for the PRTree algorithm.
    !
    ! ARGUMENTS
    !   tree [in] : tree_model, tree object.
    !
    ! Added May 2026
    !----------------------------------------------------------------------
    implicit none
    type(tree_model), intent(in) :: tree
    ! Print initialization parameters
    call labelpr("===================================== ", -1)
    call labelpr("===  Initializing PRTree model == ", -1)
    call labelpr("===================================== ", -1)
    call intpr1("    Observations (n_train):", -1, tree%n_train)
    call intpr1("    Features (n_feat):", -1, tree%n_feat)
    call intpr1("    Max terminal nodes (max_tn):", -1, tree%max_tn)
    call intpr1("    Min observations (min_obs):", -1, tree%min_obs)
    call dblepr1("    Min proportion (min_prop):", -1, tree%min_prop)
    call dblepr1("    Min probability (min_prob):", -1, tree%min_prob)
    call dblepr1("    Complexity parameter (cp):", -1, tree%cp)
    call intpr1("     Filling type (fill):", -1, tree%fill)
    call intpr1("     Proxy criterion (crit: 1 = mean, 2 = variance, 3 = both):", -1, tree%crit)
    call labelpr(" ", -1)

    ! Print distribution parameters
    call labelpr("==============================================", -1)
    call labelpr("===== Setting the distribution parameters ====", -1)
    call labelpr("==============================================", -1)
    select case (tree%dist%dist_id)
    case (1)
      call labelpr("Using the standard Gaussian distribution", -1)
    case (2)
      call dblepr1("Using the Log-normal distribution with meanlog = 0 and sdlog = ", -1, tree%dist%dist_par)
    case (3)
      call dblepr1("Using the Student's t distribution with df = ", -1, tree%dist%dist_par)
    case (4)
      call dblepr1("Using the Gamma distribution with scale = 1 and shape = ", -1, tree%dist%dist_par)
    end select
    call labelpr(" ", -1)
  end subroutine print_starting_message

  subroutine set_tree_model(tree, int_param, dble_param, keep_full, n_obs, &
    n_train, n_feat, y_train, x_train, y_test, x_test)
    !----------------------------------------------------------------------
    ! Allocates and initializes the tree structure and input data.
    ! Must be called once before building the tree.
    !
    ! ARGUMENTS
    !  tree       [inout] : tree_model, tree object to update.
    !  int_param  [in]    : Integer, parameters to build the tree
    !                       [fill, crit, max_tn, max_d, min_obs, n_cand, by_node, dist_id]
    !  dble_param [in]    : Real, parameters to build the tree
    !                       [min_prop, min_prob, cp, par]
    !  keep_full  [in]    : logical, if .true., return the full tree structure;
    !                       if .false., return only mse values.
    !  n_obs      [in]    : Integer, number of observations.
    !  n_train    [in]    : Integer, number of observations in the training data.
    !  n_feat     [in]    : Integer, number of features.
    !  y_train    [in]    : Real, response vector for training data.
    !  x_train    [in]    : Real, feature matrix for training data.
    !  y_test     [in]    : Real, response vector for testing data.
    !  x_test     [in]    : Real, feature matrix for testing data.
    !
    ! Added July/2025
    ! Last update: May/2026
    !  - added initialization for the test data
    !----------------------------------------------------------------------
    implicit none
    ! tree object
    type(tree_model), intent(inout) :: tree
    ! parameters
    integer, intent(in) :: int_param(10)
    real(dp), intent(in) :: dble_param(4)
    logical, intent(out) :: keep_full
    ! trainning and testing data
    integer, intent(in) :: n_obs, n_train, n_feat
    real(dp), intent(in), target, contiguous :: y_train(:), x_train(:, :)
    real(dp), intent(in), target, contiguous :: y_test(:)
    real(dp), intent(in), target, contiguous :: x_test(:, :)

    !--------------------------------------------------
    ! Step A: Set the tree parameters
    !--------------------------------------------------
    ! Building criteria (all)
    tree%fill = int_param(1)
    tree%crit = int_param(2)
    tree%max_tn = int_param(3)
    tree%max_d = int_param(4)
    tree%min_obs = int_param(5)
    tree%n_cand = int_param(6)
    tree%by_node = int_param(7) == 1
    tree%min_prop = dble_param(1)
    tree%min_prob = dble_param(2)
    tree%cp = dble_param(3)

    ! argsDist object
    tree%dist%dist_id = int_param(8)
    tree%dist%dist_par = dble_param(4)

    ! print info level
    tree%printinfo = int_param(9)

    ! flag for returning the full tree structure or only mse values
    keep_full = (int_param(10) == 1)

    !-----------------------------------------------------------------
    ! Step B: Set information for the training and testing data
    !-----------------------------------------------------------------
    tree%n_obs = n_obs
    tree%n_feat = n_feat
    allocate(tree%sigma(n_feat), source=0.0_dp)

    tree%n_train = n_train
    tree%y_train => y_train
    tree%x_train => x_train
    tree%na = isnan(x_train)
    tree%n_miss = count(tree%na, dim=1)
    tree%any_na = any(tree%n_miss > 0)

    tree%n_test = n_obs - n_train
    tree%x_test => x_test
    tree%y_test => y_test

    if (tree%printinfo >= log_base) call print_starting_message(tree)
  end subroutine set_tree_model

  subroutine assign_pointers(tree, keep_full, n_sigmas, yhat_train, yhat_test, &
      xregion, p_train, p_test, gammahat, nodes_info, thresholds)
    !----------------------------------------------------------------------
    ! Assigns pointers for the training and testing data, probabilities, and region assignments.
    !
    ! ARGUMENTS
    !   tree       [inout] : tree_model, tree object.
    !   keep_full  [in]    : logical, if .true., return the full tree structure;
    !                        if .false., return only mse values.
    !   n_sigmas   [in]    : Integer, number of candidates for sigma selection
    !   yhat_train [in]    : Real, predicted values for training data (n_train).
    !   yhat_test  [in]    : Real, predicted values for testing data (n_test).
    !   xregion    [in]    : Integer, assigned region for each training observation (n_train).
    !   p_train    [in]    : Real, matrix P for training data (n_train x max_tn).
    !   p_test     [in]    : Real, matrix P for testing data (n_test x max_tn).
    !   gammahat   [in]    : Real, coefficients for the tree (max_tn).
    !   nodes_info [in]    : Integer, matrix with nodes information
    !   thresholds [in]    : Real, vector of thresholds
    !
    ! Added May/2026
    !  - added pointers to matrices and vectore already allocated in R to avoid
    !    unnecessary copying and allocation in Fortran.
    !----------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout), target :: tree
    logical, intent(in) :: keep_full
    integer, intent(in) :: n_sigmas
    integer, intent(in), target, contiguous :: xregion(:)
    real(dp), intent(in), target, contiguous :: yhat_train(:)
    real(dp), intent(in), target, contiguous :: yhat_test(:)
    real(dp), intent(in), target, contiguous :: p_train(:, :)
    real(dp), intent(in), target, contiguous :: p_test(:, :)
    real(dp), intent(in), target, contiguous :: gammahat(:)
    real(dp), intent(in), target, contiguous :: thresholds(:)
    integer, intent(in), target, contiguous :: nodes_info(:, :)

    ! (a) if kepp_full = .false. R objects are dummy arguments use "temp" to build the tree.
    ! (b) if kepp_full = .true. and n_sigmas > 1 use "temp" to build and the R object to save.
    ! (c) if kepp_full = .true. and n_sigmas = 1 point directly to the R object
    if(.not. keep_full .or. n_sigmas > 1) then
      ! xregion
      allocate(tree%region_temp(tree%n_train))
      tree%region => tree%region_temp

      ! yhat_train and yhat_test
      allocate(tree%yhat_temp(tree%n_train))
      tree%yhat_train => tree%yhat_temp
      if(tree%n_test > 0) then
        allocate(tree%yhat_ttemp(tree%n_test))
        tree%yhat_test => tree%yhat_ttemp
      end if

      ! nodes_info
      allocate(tree%nodes_info_temp(2 * tree%max_tn - 1, 5))
      tree%nodes_info => tree%nodes_info_temp

      ! thresholds
      allocate(tree%thresholds_temp(2 * tree%max_tn - 1))
      tree%thresholds => tree%thresholds_temp
    else
      tree%region => xregion
      tree%yhat_train => yhat_train
      tree%yhat_test => yhat_test
      tree%nodes_info => nodes_info
      tree%thresholds => thresholds
    end if

    ! Set p_train, p_test and gammahat
    ! (a) if kepp_full = .false. R objects are dummy arguments use "temp" to build the tree.
    ! (b) if kepp_full = .true. and n_sigmas > 1 use R object as working array and "temp" as
    !     helper to save the best candidate during the search for the best sigma.
    ! (c) if kepp_full = .true. and n_sigmas = 1 point directly to the R object.
    if(.not. keep_full) then
      ! p_train and p_test
      allocate(tree%p_temp(tree%n_train, tree%max_tn))
      tree%p => tree%p_temp
      if(tree%n_test > 0) then
        allocate(tree%p_test_temp(tree%n_test, tree%max_tn))
        tree%p_test => tree%p_test_temp
      end if

      ! gammahat
      allocate(tree%gamma_temp(tree%max_tn))
      tree%gammahat => tree%gamma_temp
    else
      tree%p => p_train
      tree%p_test => p_test
      tree%gammahat => gammahat
    end if
  end subroutine assign_pointers

  !===================================================================
  ! STEP 2: Initialize the root node
  ! Initialization performed:
  !   => tree: COMPLETE  (2.1: start_root)
  !   => tree%nodes: COMPLETE (2.2: start_nodes_root)
  !   => tree%nodes%state: PARTIAL (2.3: start_node_state_root)
  !     - %imp_feature: will be initialized DURING split (n_imp = 0)
  !     - %best: will be initialized DURING split. (update = .false.)
  !===================================================================
  subroutine start_nodes_root(node, n_train, n_feat, max_tn)
    !-------------------------------------------------------------
    !  Allocates and initializes the root node and its bounds.
    !
    ! ARGUMENTS
    !   node   [inout] : tree_node(:), array of nodes.
    !   n_train[in]    : Integer, number of observations.
    !   n_feat [in]    : Integer, number of features.
    !   max_tn [in]    : Integer, maximum number of terminal nodes.
    !
    ! Details:
    ! Some variables will require initialization/update before use
    !   - imp_feature: will be initialized DURING split (n_imp = 0)
    !   - thresholds: will be updated BEFORE split
    !   - best: will be initialized DURING split. (update = .false.)
    !
    ! Added July/2025
    !-------------------------------------------------------------
    implicit none
    type(tree_node), allocatable, intent(out) :: node(:)
    integer, intent(in) :: n_train, n_feat, max_tn
    integer :: i

    ! Allocate the nodes object and the bounds argument
    allocate(node(2 * max_tn - 1))

    ! initializes the state variables
    node(1)%split = .true.
    node(1)%update = .true.    ! need update
    node(1)%n_cand_found = 0   ! no candidates found at this point
    node(1)%n_obs_node = n_train ! all data starts in the root node
    node(1)%n_imp = 0          ! to start the update
    node(1)%idx = [(i, i=1, n_train)]
    allocate(node(1)%imp_feat(n_feat))
    node(1)%imp_feat = 0       ! to start the update
    allocate(node(1)%thresholds(n_feat))
  end subroutine start_nodes_root

  subroutine start_root(tree)
    !--------------------------------------------------------------
    ! Initializes the tree structure and root node.
    ! Initialization performed:
    !   => tree: COMPLETE  (2.1: start_root)
    !   => tree%nodes: COMPLETE (2.2: start_nodes_root)
    !   => tree%nodes%state: PARTIAL (2.3: start_node_state_root)
    !     - %imp_feature: will be initialized DURING split (n_imp = 0)
    !     - %best: will be initialized DURING split. (update = .false.)
    !
    ! ARGUMENTS
    !   tree     [inout] : tree_model, tree object.
    !
    ! Added July/2025
    !--------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer :: f, n
    logical, allocatable :: x_mask(:)

    ! At start number of nodes = number of terminal nodes
    tree%n_nodes = 1
    tree%n_tn = 1
    tree%n_cand_found = 0

    ! In the root node
    ! - Region = 1 and tn_id = 1
    ! - P = 1 and gammahat = mean(y)
    ! Use vector so there is no need to check if the object is already allocated
    tree%region = 1
    if (allocated(tree%tn_id)) deallocate(tree%tn_id)
    allocate(tree%tn_id(tree%max_tn))
    tree%tn_id(1) = 1
    tree%p = 0.0_dp
    tree%p(:, 1) = 1.0_dp
    tree%gammahat = 0.0_dp
    tree%gammahat(1) = sum(tree%y_train) / tree%n_train
    tree%yhat_train = tree%gammahat(1)
    tree%mse_train = sum(tree%y_train**2) - tree%n_train * tree%gammahat(1)**2

    ! initialize the nodes_info matrix for the root
    tree%nodes_info = 0
    tree%nodes_info(1, idx_ni_id) = 1
    tree%nodes_info(1, idx_ni_terminal) = 1
    tree%nodes_info(1, idx_ni_father) = 0
    tree%nodes_info(1, idx_ni_depth) = 0
    tree%nodes_info(1, idx_ni_feature) = 0

    ! initialize the thresholds array
    tree%thresholds = 0.0_dp

    ! Allocate the root node and its components
    call start_nodes_root(tree%nodes, tree%n_train, tree%n_feat, tree%max_tn)

    ! Loop to compute the thresholds for each feature
    !  - checks if the node can be splitted using min_obs
    if (tree%printinfo >= log_debug) then
      ! Debug: thresholds
      call labelpr("------------------------------------------------------------------", -1)
      call labelpr("    Getting the thresholds for the features at the root node", -1)
      call labelpr("------------------------------------------------------------------", -1)
    end if

    allocate(x_mask(tree%n_train))
    loop_find_thr: do f = 1, tree%n_feat
      ! Create mask for non-missing values of feature f
      if (tree%n_miss(f) == 0) then
        x_mask = .true.
      else
        x_mask = .not. tree%na(:, f)
      end if
      n = tree%n_train - tree%n_miss(f)

      ! Skip feature if insufficient non-missing values
      if (n < 2 * tree%min_obs) then
        tree%nodes(1)%thresholds(f)%nt = 0
        cycle
      end if

      ! Find thresholds for current feature (considering only non-missing values)
      call find_thresholds(tree%min_obs, tree%x_train, x_mask, f, tree%nodes(1)%thresholds(f))

      if (tree%printinfo >= log_debug) then
        ! Debug: thresholds found
        call intpr1("    Feature ", -1, f)
        if (tree%nodes(1)%thresholds(f)%nt == 0) then
          call labelpr("    No thresholds found", -1)
        else
          call intpr1("    Number of thresholds:", -1, tree%nodes(1)%thresholds(f)%nt)
        end if
        call labelpr(" ", -1)
      end if
    end do loop_find_thr
  end subroutine start_root

  !===============================================================================
  ! STEP 3: Build the three
  ! Loop to grow the three
  !   3.1 update_split_candidates:
  !       Set/Update split candidates (depth, min_prop, min_prob and max_d)
  !      (last two child nodes plus the ones not discarded during split)
  !   3.2 update_state_thresholds:
  !      If tn > 1, update thresholds for child nodes
  !      If tn > 1 update the important features for child nodes
  !   3.3 find_node_splits:
  !       Find the best n_candidates by node.
  !       (each node only goes through the process once)
  !   3.4 split_full_analysis:
  !       Perform the full analysis using the best candidates.
  !       (need to be done every loop because it depends on P)
  !===============================================================================
  pure subroutine update_node_to_terminal(node)
    !-------------------------------------------------------------
    ! Marks a node as terminal and disables further splits.
    !
    ! ARGUMENTS
    !   node [inout] : tree_node, node to mark as terminal.
    !
    ! Added July/2025
    !-------------------------------------------------------------
    implicit none
    type(tree_node), intent(inout) :: node

    ! disable further splits
    node%split = .false.

    ! Information on splitting is no longer required.
    ! Deallocate the components from state variable
    node%n_cand_found = 0
    node%n_imp = 0
    if (allocated(node%idx)) deallocate(node%idx)
    if (allocated(node%imp_feat)) deallocate(node%imp_feat)
    if (allocated(node%thresholds)) deallocate(node%thresholds)
    if (allocated(node%best)) deallocate(node%best)
  end subroutine update_node_to_terminal

  pure subroutine update_split_candidates(tree, n_split, do_split)
    !---------------------------------------------------------------------------
    ! Checks which nodes are eligible for splitting based on the depth of
    ! the node and probability criteria. Number of observations are checked
    ! in another step
    !
    ! ARGUMENTS
    !   tree     [inout] : tree_model, tree structure.
    !   n_split  [inout] : Integer, number of splittable nodes (old vs current).
    !   do_split [inout] : Integer, indexes of splittable nodes (old vs current).
    !
    ! Added July/2025
    ! Last update: May/2026
    !   - now passes tree instead of net, ctrl and n_obs
    !---------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(inout) :: n_split
    integer, allocatable, intent(inout) :: do_split(:)
    ! local variables
    real(dp) :: perc_comp, multi, min_prob
    integer :: j, n_found, tn_id, nodeid
    integer, allocatable :: cut(:)

    ! copy the id of the previous splittable nodes
    allocate(cut(n_split + 2))
    cut = [do_split, tree%n_nodes - 1, tree%n_nodes]

    ! check the current split condition of the node
    !  - split condition can change during the search for thresholds
    !    when the number of observations is evaluated.
    !  - the node is updated to terminal when the split condition changes.

    multi = 1.0_dp / real(tree%n_train, dp)
    min_prob = 2.0_dp * tree%min_prob
    n_found = 0
    do j = 1, n_split + 2
      nodeid = cut(j)

      if (.not. tree%nodes(nodeid)%split) cycle ! node updated to terminal during the process

      ! Check stopping criteria (only for last two child nodes):
      !  - depth and the percentage of probabilities higher than a threshold
      !  - use P > 2 * min_prob as a proxy for the next division.
      if (j > n_split) then
        if (tree%nodes_info(nodeid, idx_ni_depth) >= tree%max_d) then
          call update_node_to_terminal(tree%nodes(nodeid))
          cycle
        end if

        ! P(child) < a, whenever P(father) < a.
        ! If a <= P(father) < 2*a implies P(child1) < a, whenever P(child2) > a
        tn_id = merge(tree%n_tn - 1, tree%n_tn, j == n_split + 1)
        perc_comp = count(tree%p(:, tn_id) > min_prob) * multi
        if (perc_comp <= tree%min_prop) then
          call update_node_to_terminal(tree%nodes(nodeid))
          cycle
        end if
      end if

      n_found = n_found + 1 ! found a candidate to split
      cut(n_found) = nodeid ! id of split candidate
    end do

    ! save the position of cuttable nodes
    n_split = n_found
    if (n_split > 0) do_split = cut(1:n_split)
  end subroutine update_split_candidates

  subroutine print_message(left, f, n, printinfo)
    !----------------------------------------------------------------------
    ! Helper subroutine to print debug messages about thresholds found
    ! for a given feature in a node.
    !
    ! ARGUMENTS
    !   left [in] : logical, .true. if left node, .false. if right node.
    !   f    [in] : Integer, feature index.
    !   n    [in] : Integer, number of thresholds found.
    !   printinfo [in] : Integer, log level.
    !
    ! Added July/2025
    !----------------------------------------------------------------------
    implicit none
    logical, intent(in) :: left
    integer, intent(in) :: f, n, printinfo
    if (printinfo >= log_debug) then
      ! Debug: thresholds found
      if (left) then
        call intpr1(" ---- Thresholds for splitting feature (left node) ----", -1, f)
      else
        call intpr1(" ---- Thresholds for splitting feature (right node) ----", -1, f)
      end if
      if (n == 0) then
        call labelpr("    No thresholds found for this feature", -1)
      else
        call intpr1("    Number of thresholds found for this feature:", -1, n)
      end if
      call labelpr(" ", -1)
    end if
  end subroutine print_message

  subroutine update_state_thresholds(tree, updated, fa_id, l_id, r_id)
    !----------------------------------------------------------------------
    ! Updates thresholds for child nodes after a split operation.
    ! Splits the threshold vector for current variable and recompute for others
    !
    ! ARGUMENTS
    !   tree    [inout] : tree_model, the tree structure.
    !   updated [inout] : logical(2), if .true., new splits can be attempted
    !   fa_id   [in]    : Integer, id of the father node.
    !   l_id    [in]    : Integer, id of the left child node.
    !   r_id    [in]    : Integer, id of the right child node.
    !
    ! NOTES
    !   Thresholds are allocated during split and only updated here
    !
    !   For current variable:
    !    - Checks min_obs requirement for future splits in both directions
    !    - Preserves original threshold ordering
    !
    ! Added July/2025
    ! Last update: May/2026
    !  removed the arguments from the parent node and replaced with the tree
    !  structure to avoid slicing x and na.
    !----------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    logical, intent(inout) :: updated(2)
    integer, intent(in) :: fa_id, l_id, r_id

    ! Local variables
    integer :: pos, i, f, feat, n_node, n_obs_node
    real(dp) :: thr
    logical, allocatable :: msk_left(:)
    logical, allocatable :: msk_right(:)
    logical, allocatable :: mask_miss(:)
    type(info_thr) :: feat_thr
    integer :: n1, n2, n, n_train, ii
    logical, allocatable :: x_mask(:)

    n_train = tree%n_train
    n_obs_node = tree%nodes(fa_id)%n_obs_node

    ! masks for the observations in the left and right child nodes
    ! intially set to .false. and updated to .true. for the observations in the node
    !  - idx gives the indexes that need update
    !  - best(1)%regid gives the region of the left child node (after split)
    allocate(msk_left(n_train), source = .false.)
    allocate(msk_right(n_train), source = .false.)
    do i = 1, n_obs_node
      ii = tree%nodes(fa_id)%idx(i)
      msk_left(ii) = tree%nodes(fa_id)%best(1)%regid(i) == left_id
      msk_right(ii) = .not. msk_left(ii)
    end do

    allocate(mask_miss(n_train)) ! mask for non-missing values in the node
    allocate(x_mask(n_train))    ! mask for x values in the node (non-missing and in the l/r child node)

    if (tree%printinfo >= log_debug) then
      ! Debug: thresholds found
      call labelpr("---------------------------------------------", -1)
      call labelpr("    Updating thresholds for the new nodes", -1)
      call labelpr("---------------------------------------------", -1)
    end if

    ! splitting feature and threshold
    feat = tree%nodes_info(fa_id, idx_ni_feature)
    thr = tree%thresholds(fa_id)

    ! counts the number of thresholds found globally (left and right)
    n1 = 0
    n2 = 0

    if (tree%nodes(l_id)%split) tree%nodes(l_id)%thresholds%nt = 0
    if (tree%nodes(r_id)%split) tree%nodes(r_id)%thresholds%nt = 0

    do f = 1, tree%n_feat
      ! current feature will be processed at the end of the loop
      if (f == feat) cycle

      ! mask for non-missing values in the current node
      mask_miss =  .not. tree%na(:, f)

      ! recompute thresholds for left node
      n = 0         ! local counter
      if (tree%nodes(l_id)%split) then
        x_mask = msk_left .and. mask_miss
        call find_thresholds(tree%min_obs, tree%x_train, x_mask, f, tree%nodes(l_id)%thresholds(f))
        n = tree%nodes(l_id)%thresholds(f)%nt
        n1 = n1 + n
      end if
      call print_message(.true., f, n, tree%printinfo)

      ! recompute thresholds for right node
      n = 0
      if (tree%nodes(r_id)%split) then
        x_mask = msk_right .and. mask_miss
        call find_thresholds(tree%min_obs, tree%x_train, x_mask, f, tree%nodes(r_id)%thresholds(f))
        n = tree%nodes(r_id)%thresholds(f)%nt
        n2 = n2 + n
      end if
      call print_message(.false., f, n, tree%printinfo)
    end do

    ! Finds the position of thr in the thresholds vector
    pos = findloc(tree%nodes(fa_id)%thresholds(feat)%thr >= thr, value=.true., dim=1)

    ! mask for non missing values
    mask_miss =  .not. tree%na(:, feat)

    ! thresholds information for splitting feature
    feat_thr = tree%nodes(fa_id)%thresholds(feat)

    ! LEFT NODE: can only split using thresholds < thr
    ! Since these observations are already in left node (msk = .true.),
    ! we only need to ensure the RIGHT child after split has >= min_obs
    n = 0   ! local counter
    if (tree%nodes(l_id)%split .and. pos > 1) then ! There exist at least one thresholds < thr
      do i = pos - 1, 1, -1                 ! Check from largest < thr downward
        n_node = count(tree%x_train(tree%nodes(l_id)%idx, feat) > feat_thr%thr(i) .and. &
                       .not. tree%na(tree%nodes(l_id)%idx, feat))
        if (n_node >= tree%min_obs) then
          tree%nodes(l_id)%thresholds(feat)%nt = i
          tree%nodes(l_id)%thresholds(feat)%thr = feat_thr%thr(1:i)
          exit
        end if
      end do
      n = tree%nodes(l_id)%thresholds(feat)%nt
      n1 = n1 + n
    end if
    call print_message(.true., feat, n, tree%printinfo)

    ! RIGHT NODE: can only split using thresholds > thr
    ! Since these observations are already in right node (msk = .false.),
    ! we only need to ensure the LEFT child after split has >= min_obs
    n = 0  ! local counter
    if (tree%nodes(r_id)%split .and. pos < feat_thr%nt) then ! There exist at least one thresholds > thr
      do i = pos + 1, feat_thr%nt                      ! Start after the exact split threshold
        n_node = count(tree%x_train(tree%nodes(r_id)%idx, feat) <= feat_thr%thr(i) .and. &
                       .not. tree%na(tree%nodes(r_id)%idx, feat))
        if (n_node >= tree%min_obs) then
          tree%nodes(r_id)%thresholds(feat)%nt = feat_thr%nt - i + 1
          tree%nodes(r_id)%thresholds(feat)%thr = feat_thr%thr(i:feat_thr%nt)
          exit
        end if
      end do
      n = tree%nodes(r_id)%thresholds(feat)%nt
      n2 = n2 + n
    end if
    call print_message(.false., feat, n, tree%printinfo)

    ! check if the nodes can be splitted
    ! check if the nodes can be splitted
    updated(1) = (n1 > 0)
    updated(2) = (n2 > 0)
  end subroutine update_state_thresholds

  subroutine update_node_state(tree, fail)
    !----------------------------------------------------------------------
    ! Updates the node state for child nodes after a split.
    !
    ! ARGUMENTS
    !   tree [inout] : tree_model, the tree structure.
    !   fail [inout] : logical, if .false. no threholds were found.
    !
    ! Added July/2025
    ! Last update: May/2026
    !   removed the arguments from the parent node and replaced with the tree structure
    !   since the parent node is always the last splitted node, we can access it through
    !   the tree structure. This also avoid slicing x and na
    !----------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    logical, intent(inout) :: fail
    logical :: update(2)
    integer :: fa_id, l_id, r_id
    integer :: n_nodes
    integer :: i, c_l, c_r

    ! number of nodes for the current tree structure (after split)
    n_nodes = tree%n_nodes

    ! Get the id of the father node: the last two nodes are the child nodes
    ! of the last splitted node
    fa_id = tree%nodes_info(n_nodes, idx_ni_father)
    l_id = n_nodes - 1
    r_id = n_nodes

    if (.not. (tree%nodes(l_id)%split .or. tree%nodes(r_id)%split)) then
      call update_node_to_terminal(tree%nodes(l_id))
      call update_node_to_terminal(tree%nodes(r_id))
      call update_node_to_terminal(tree%nodes(fa_id))
      fail = .true.
      return
    end if

    if (tree%nodes(l_id)%split) then
      tree%nodes(l_id)%n_obs_node = tree%nodes(fa_id)%best(1)%n_l  ! size
      allocate(tree%nodes(l_id)%thresholds(tree%n_feat))           ! thresholds
      allocate(tree%nodes(l_id)%imp_feat(tree%n_feat), source = 0) ! pre-allocate to max capacity
    end if

    if (tree%nodes(r_id)%split) then
      tree%nodes(r_id)%n_obs_node = tree%nodes(fa_id)%best(1)%n_r  ! size
      allocate(tree%nodes(r_id)%thresholds(tree%n_feat))           ! thresholds
      allocate(tree%nodes(r_id)%imp_feat(tree%n_feat), source = 0) ! pre-allocate to max capacity
    end if

    if (tree%nodes(l_id)%split .and. tree%nodes(r_id)%split) then
      allocate(tree%nodes(l_id)%idx(tree%nodes(l_id)%n_obs_node))
      allocate(tree%nodes(r_id)%idx(tree%nodes(r_id)%n_obs_node))
      c_l = 1
      c_r = 1
      do i = 1, tree%nodes(fa_id)%n_obs_node
        if (tree%nodes(fa_id)%best(1)%regid(i) == left_id) then
          tree%nodes(l_id)%idx(c_l) = tree%nodes(fa_id)%idx(i)
          c_l = c_l + 1
        else
          tree%nodes(r_id)%idx(c_r) = tree%nodes(fa_id)%idx(i)
          c_r = c_r + 1
        end if
      end do
    else if (tree%nodes(l_id)%split) then
      tree%nodes(l_id)%idx = pack(tree%nodes(fa_id)%idx, mask=tree%nodes(fa_id)%best(1)%regid == left_id)
    else if (tree%nodes(r_id)%split) then
      tree%nodes(r_id)%idx = pack(tree%nodes(fa_id)%idx, mask=tree%nodes(fa_id)%best(1)%regid == right_id)
    end if

    call update_state_thresholds(tree, update, fa_id, l_id, r_id)

    ! If no thresholds were found, deallocate the old node state to free memory
    fail = .not. (update(1) .or. update(2))
    if (fail) then
      call update_node_to_terminal(tree%nodes(fa_id))
      call update_node_to_terminal(tree%nodes(l_id))
      call update_node_to_terminal(tree%nodes(r_id))
      return
    end if

    ! At least one child has valid thresholds, update the NA masks.
    ! Copy the important features from the parent node (if needed)
    if (update(1)) then
      tree%nodes(l_id)%n_imp = tree%nodes(fa_id)%best(1)%n_imp
      tree%nodes(l_id)%imp_feat = tree%nodes(fa_id)%best(1)%imp_feat
    else
      call update_node_to_terminal(tree%nodes(l_id))
    end if

    if (update(2)) then
      tree%nodes(r_id)%n_imp = tree%nodes(fa_id)%best(1)%n_imp
      tree%nodes(r_id)%imp_feat = tree%nodes(fa_id)%best(1)%imp_feat
    else
      call update_node_to_terminal(tree%nodes(r_id))
    end if

    ! Reset the old node state to free memory
    call update_node_to_terminal(tree%nodes(fa_id))
  end subroutine update_node_state

  pure subroutine update_best_list(n_cand, best, try, fail)
    !-----------------------------------------------------------------------------------
    !  Updates the list of best candidate splits for a node.
    !
    ! ARGUMENTS
    !   n_cand [in]    : Integer, number of candidates.
    !   best   [inout] : info_split, current best candidates.
    !   try    [in]    : info_split, new candidate.
    !   fail   [inout] : Logical, flag indicating if update failed.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    integer, intent(in) :: n_cand
    type(info_split), intent(inout) :: best(:)
    type(info_split), intent(in) :: try
    logical, intent(inout) :: fail
    integer :: worst_pos

    fail = .true.
    worst_pos = minloc(best(1:n_cand)%score, dim=1)
    if (best(worst_pos)%score < try%score) then
      best(worst_pos) = try
      fail = .false.
    end if
  end subroutine update_best_list

  pure subroutine update_best_list_global(n_best, best_id, best_score, nodeid, &
    n_try, try, by_node, first, n_cand_found)
    !-----------------------------------------------------------------------------------
    ! Updates the global list of best candidate splits across all nodes.
    !
    ! ARGUMENTS
    !   n_best       [in]    : Integer, maximum number of best candidates.
    !   best_id      [inout] : integer, current best candidates.
    !   best_score   [inout] : real, current best scores.
    !   nodeId       [in]    : Integer, node id.
    !   n_try        [in]    : Integer, number of new candidates.
    !   try          [in]    : info_split, new candidates.
    !   by_node      [in]    : Logical, flag for node-wise update.
    !   first        [inout] : Logical, flag for first update.
    !   n_cand_found [inout] : Integer, number of candidates found.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    integer, intent(in) :: n_try, nodeid, n_best
    integer, intent(inout) :: best_id(:, :)
    real(dp), intent(inout) :: best_score(:)
    type(info_split), intent(in) :: try(:)
    logical, intent(in) :: by_node
    logical, intent(inout) :: first
    integer, intent(inout) :: n_cand_found
    integer :: i, worst_pos
    logical :: isfirst

    isfirst = first
    if (isfirst) then
      n_cand_found = 0
      best_id = 0
      best_score = 0.0_dp
      first = .false.
    end if

    ! if by_node, merge the best candidates
    if (isfirst .or. by_node) then
      ! merge the best candidates with the new ones. No need to sort.
      best_id(n_cand_found + 1:n_cand_found + n_try, 1) = nodeid
      best_id(n_cand_found + 1:n_cand_found + n_try, 2) = [(i, i=1, n_try)]
      best_score(n_cand_found + 1:n_cand_found + n_try) = try%score
      n_cand_found = n_cand_found + n_try
      return
    end if

    ! if by_node = .false., find the worst candidate and update it
    ! with the new candidate if it is better
    do i = 1, n_try
      ! find the worst candidate position
      worst_pos = minloc(best_score(1:n_best), dim=1)
      if (try(i)%score > best_score(worst_pos)) then
        best_score(worst_pos) = try(i)%score
        best_id(worst_pos, 1) = nodeid
        best_id(worst_pos, 2) = i
        n_cand_found = n_cand_found + 1
      end if
    end do
  end subroutine update_best_list_global

  pure subroutine get_idx_prob(idx_p, na, n_imp, imp_feat, n_prob, idx, fill, f)
    !---------------------------------------------------------------------------------
    ! Updates the index information for probability computation in the current node.
    ! This subroutine determines which observations should be used to compute the
    ! probability matrix `P` for the left and right child nodes, based on the
    ! missingness of features and the chosen fill strategy.
    !
    ! ARGUMENTS
    !   idx_p    [inout] : idx_prob, information on missing values.
    !   n_train  [in]    : Integer, total number of observations.
    !   na       [in]    : Logical, complete missing mask matrix.
    !   n_imp    [in]    : Integer, number of important features.
    !   imp_feat [in]    : Integer, indices of important features.
    !   n_prob   [in]    : Integer, number of observations with probabilities > eps.
    !   idx      [in]    : Integer, global indices with p > eps.
    !   fill     [in]    : Integer, fill strategy for missing values.
    !   f        [in]    : Integer, splitting feature
    !
    ! Added: July/2025
    ! Last update: March/2026
    !   - replaced n_miss(2) with n_miss_f and n_miss_any
    !---------------------------------------------------------------------------------
    implicit none
    type(idx_prob), intent(inout) :: idx_p
    integer, intent(in) :: n_imp, n_prob, fill, f
    logical, intent(in) :: na(:, :)
    integer, intent(in) :: imp_feat(:)
    integer, intent(in) :: idx(:)
    logical, allocatable :: any_na(:)
    integer :: j

    idx_p%n_prob = n_prob

    if(fill < 2) then
      allocate(any_na(n_prob), source = .false.)
      do j = 1, n_imp
        any_na = any_na .or. na(idx, imp_feat(j))
      end do
    end if

    ! Set the indices to be used to update the probability matrix P
    select case (fill)
    case (0)
      ! Only depends on whether the observation has any missing values or not
      idx_p%n_miss_any = count(any_na)  ! any missing
      idx_p%n_miss_f = 0  ! dummy value
      if (idx_p%n_miss_any > 0) then
        idx_p%idx_any_na = pack(idx, any_na)
      end if
      if (idx_p%n_miss_any < n_prob) then
        idx_p%idx = pack(idx,.not. any_na)
      end if
    case (1)
      ! Depends on both: whether the observation has any missing values or not
      ! and on whether the splitting feature is missing
      idx_p%n_miss_any = count(any_na)  ! any missing
      idx_p%n_miss_f = count(na(idx, f)) ! Xf missing
      if (idx_p%n_miss_f > 0) then
        idx_p%idx_na_f = pack(idx, na(idx, f))
      end if
      if (idx_p%n_miss_f < idx_p%n_miss_any) then
        idx_p%idx_any_na = pack(idx, any_na .and. (.not. na(idx, f)))
      end if
      if (idx_p%n_miss_any < n_prob) then
        idx_p%idx = pack(idx,.not. any_na)
      end if
    case (2)
      ! Only depends on whether the splitting feature is missing or not
      idx_p%n_miss_any = 0  ! dummy value
      idx_p%n_miss_f = count(na(idx, f))
      if (idx_p%n_miss_f > 0) then
        idx_p%idx_na_f = pack(idx, na(idx, f))
      end if
      if (idx_p%n_miss_f < n_prob) then
        idx_p%idx = pack(idx,.not. na(idx, f))
      end if
    end select
  end subroutine get_idx_prob

  pure subroutine assign_missing(info, n_miss, y, y_sum, idx, crit)
    !---------------------------------------------------------------------------
    ! Assigns a set of observations with missing feature values to left or
    ! right node using different criteria based on input selection.
    !
    ! ARGUMENTS
    !   info    [inout] : info_split, split information.
    !   n_miss  [in]    : Integer, number of observations (local).
    !   y       [in]    : Real, response values for missing (local).
    !   y_sum   [inout] : Real, sum for left and right nodes (local).
    !   idx     [in]    : Integer, index of the missing feature (local).
    !   crit    [in]    : Integer, selection of assignment criterion:
    !                      1 - Maximize difference in means
    !                      2 - Maximize between-node variability
    !                      3 - Sum of both criteria
    !
    ! Added July/2025
    !---------------------------------------------------------------------------
    implicit none
    type(info_split), intent(inout) :: info
    integer, intent(in) :: n_miss
    real(dp), intent(in) :: y(:)
    real(dp), intent(inout) :: y_sum(2)
    integer, intent(in) :: idx(:), crit
    real(dp) :: score(2)
    integer :: new_nl, new_nr
    real(dp) :: new_sl, new_sr
    integer :: i

    do i = 1, n_miss
      ! Pre-calculate updated sums and counts to avoid redundant calculations
      new_sl = y_sum(1) + y(i)
      new_nl = info%n_l + 1
      new_sr = y_sum(2) + y(i)
      new_nr = info%n_r + 1
      score = 0.0_dp

      ! Calculate score based on the selected criterion
      if (crit == 1 .or. crit == 3) then
        ! Maximize difference in means
        score(1) = abs((new_sl / new_nl) - (y_sum(2) / info%n_r))
        score(2) = abs((y_sum(1) / info%n_l) - (new_sr / new_nr))
      end if
      if (crit == 2 .or. crit == 3) then
        ! Maximize between-node variability (sum of squares)
        score(1) = score(1) + new_sl**2 / new_nl + y_sum(2)**2 / info%n_r
        score(2) = score(2) + y_sum(1)**2 / info%n_l + new_sr**2 / new_nr
      end if

      ! Assign to maximize the score
      if (score(1) > score(2)) then
        y_sum(1) = new_sl
        info%n_l = new_nl
        info%regid(idx(i)) = left_id
      else
        y_sum(2) = new_sr
        info%n_r = new_nr
        info%regid(idx(i)) = right_id
      end if
    end do
  end subroutine assign_missing

  subroutine update_prob_no_na(argsd, n_idx, idx, x_f, threshold, &
    b_lower, b_upper, sigma_f, prob, p)
    implicit none
    !----------------------------------------------------------------------------
    ! Computes the probability for the left and right node given a feature and a
    ! threshold when there are no missing values in the important features.
    ! WARNING: This subrotuine implicitly assumes that for idx, all prob > eps.
    !
    ! ARGUMENTS:
    !   argsd     [in] : argsDist, distribution related parameters.
    !   n_idx     [in] : Integer, size of idx.
    !   idx       [in] : Integer, global indices to compute P(idx, :).
    !   x_f       [in] : Real, the value for the splitting feature.
    !   threshold [in] : Real, the threshold for the splitting feature.
    !   b_lower   [in] : Real, lower bound for the splitting feature.
    !   b_upper   [in] : Real, upper bound for the splitting feature.
    !   sigma_f   [in] : Real, sigma for splitting feature.
    !   prob      [in] : probability for the father node.
    !   p        [out] : probability for the left and right nodes.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n_idx
    integer, intent(in) :: idx(:)
    real(dp), intent(in) :: x_f(:)
    real(dp), intent(in) :: threshold, b_lower, b_upper, sigma_f
    real(dp), intent(in) :: prob(:)
    real(dp), intent(inout) :: p(:, :)
    real(dp), allocatable :: p_temp(:)
    real(dp) :: p_sum
    integer :: i

    allocate(p_temp(n_idx))
    ! Compute P([-Inf, threshold])
    p_temp = pdist(argsd, n_idx, threshold, x_f(idx), sigma_f)

    ! If the lower bound is -Inf, the left node covers (-Inf, threshold]
    ! Otherwise, the left node covers [lower, threshold]
    if (b_lower <= neg_inf) then
      p(idx, 1) = p_temp
    else
      p(idx, 1) = p_temp - pdist(argsd, n_idx, b_lower, x_f(idx), sigma_f)
    end if

    ! If the upper bound is Inf, the right node covers (threshold, Inf)
    ! Otherwise, the right node covers [threshold, upper]
    if (b_upper >= pos_inf) then
      p(idx, 2) = 1.0_dp - p_temp
    else
      p(idx, 2) = pdist(argsd, n_idx, b_upper, x_f(idx), sigma_f) - p_temp
    end if

    ! Standardize the probabilities so that they sum to P(Parent)
    do i = 1, n_idx
      p_sum = p(idx(i), 1) + p(idx(i), 2)
      if (p_sum > eps) then
        p_temp(i) = prob(idx(i)) / p_sum
        p(idx(i), 1) = p(idx(i), 1) * p_temp(i)
        p(idx(i), 2) = p(idx(i), 2) * p_temp(i)
      else
        ! Fallback: assign deterministically based on threshold
        if (x_f(idx(i)) <= threshold) then
          p(idx(i), 1) = prob(idx(i))
          p(idx(i), 2) = 0.0_dp
        else
          p(idx(i), 1) = 0.0_dp
          p(idx(i), 2) = prob(idx(i))
        end if
      end if
    end do
  end subroutine update_prob_no_na

  subroutine update_prob(argsd, idx_p, x_f, threshold, b_lower, b_upper, sigma_f, fill, prob, p)
    !---------------------------------------------------------------------------------------
    ! Computes the probability for the left and right node given a feature and a threshold.
    ! WARNING: This subroutine implicitly assumes that for idx, all prob > eps.
    !
    ! ARGUMENTS
    !   argsd    [in] : argsDist, distribution related parameters.
    !   idx_p    [in] : type(idx_prob), information on missing values.
    !   x_f      [in] : Real, the value for the splitting feature.
    !   na_f     [in] : Logical, missing mask for the current feature.
    !   bounds   [in] : Real, lower and upper bounds for the splitting feature.
    !   sigma_f  [in] : Real, sigma for splitting feature.
    !   fill     [in] : Integer, fill strategy for missing values.
    !   n_miss   [in] : Integer, number of rows with missing values and P > eps
    !                   (global and for feature f)
    !   prob     [in] : probability for the father node.
    !   p       [out] : probability for the left and right nodes.
    !
    ! DETAILS
    !
    ! The probability is computed as
    !   P(Left)  = P(R_f_left | X_f, sigma_f) * P(R_{-f} | X_{-f}, sigma_{-f})
    !   P(Right) = P(R_f_right | X_f, sigma_f) * P(R_{-f} | X_{-f}, sigma_{-f})
    ! where
    !   P(R_{-f} | X_{-f}, sigma_{-f}) = prod_{j != f} P(R_j | X_j, sigma_j)
    !
    ! Moreover, for the parent node we have
    ! P(Parent) = prod_{j} P(R_j | X_j, sigma_j)
    !           = P(R_f | X_f, sigma_f) * prod_{j != f} P(R_j | X_j, sigma_j)
    ! Thus,
    !   P(R_{-f} | X_{-f}, sigma_{-f}) = P(Parent) / P(R_f | X_f, sigma_f)
    ! and we can rewrite
    !   P(Left)  = P(R_f_left | X_f, sigma_f) * P(Parent) / P(R_f | X_f, sigma_f)
    !   P(Right) = P(R_f_right | X_f, sigma_f) * P(Parent) / P(R_f | X_f, sigma_f)
    ! with
    !   P(R_f | X_f, sigma_f) = P(R_f_left | X_f, sigma_f) + P(R_f_right | X_f, sigma_f)
    !
    ! This means we only need to compute the marginal probabilities for the splitting feature
    ! and then scale them by the probability of the parent node.
    !
    ! Added July/2025
    !---------------------------------------------------------------------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    type(idx_prob), intent(in) :: idx_p
    integer, intent(in) :: fill
    real(dp), intent(in) :: x_f(:), sigma_f, threshold
    real(dp), intent(in) :: b_lower, b_upper
    real(dp), intent(in) :: prob(:)
    real(dp), intent(inout) :: p(:, :)
    integer :: i, n_na

    select case (fill)
    case (0)
      ! fill = 0: Assing uniform probability for both nodes if any observation in Xi is missing
      if (idx_p%n_miss_any > 0) then
        p(idx_p%idx_any_na, 1) = prob(idx_p%idx_any_na) * 0.5_dp
        p(idx_p%idx_any_na, 2) = prob(idx_p%idx_any_na) * 0.5_dp
        if (idx_p%n_miss_any == idx_p%n_prob) return
      end if
    case (1)
      ! fill = 1: Assigns uniform probability for both nodes if x_f is missing
      !           Assign 0/1 weights when X_f is not missing
      ! CASE 1: Missing values for X_f
      if (idx_p%n_miss_f > 0) then
        p(idx_p%idx_na_f, 1) = prob(idx_p%idx_na_f) * 0.5_dp
        p(idx_p%idx_na_f, 2) = prob(idx_p%idx_na_f) * 0.5_dp
        if (idx_p%n_miss_f == idx_p%n_prob) return
      end if
      ! CASE 2: X_f observed but missing values in X
      if (idx_p%n_miss_f < idx_p%n_miss_any) then
        n_na = idx_p%n_miss_any - idx_p%n_miss_f
        do i = 1, n_na
          if (x_f(idx_p%idx_any_na(i)) > threshold) then
            ! probability in the right node
            p(idx_p%idx_any_na(i), 1) = 0.0_dp
            p(idx_p%idx_any_na(i), 2) = prob(idx_p%idx_any_na(i))
          else
            ! probability in the right node
            p(idx_p%idx_any_na(i), 1) = prob(idx_p%idx_any_na(i))
            p(idx_p%idx_any_na(i), 2) = 0.0_dp
          end if
        end do
        if (idx_p%n_miss_any == idx_p%n_prob) return
      end if
    case (2)
      ! fill = 2: Assigns probability based weights.
      ! Compute the indexes for missing values using na_f
      if (idx_p%n_miss_f > 0) then
        p(idx_p%idx_na_f, 1) = prob(idx_p%idx_na_f) * 0.5_dp
        p(idx_p%idx_na_f, 2) = prob(idx_p%idx_na_f) * 0.5_dp
        if (idx_p%n_miss_f == idx_p%n_prob) return
      end if
    end select

    ! Compute the probability for complete cases (also for fill = 2 and x_f not missing)
    call update_prob_no_na(argsd, size(idx_p%idx), idx_p%idx, x_f, &
      threshold, b_lower, b_upper, sigma_f, prob, p)
  end subroutine update_prob

  subroutine update_pmatrix(tree, nodeId, splitId, fa_id, p_out)
    !---------------------------------------------------------------------------------------
    ! Updates the probability matrix for the child nodes after a split.
    !
    ! ARGUMENTS
    !   tree     [in] : tree_model, the tree structure.
    !   nodeId   [in] : Integer, node id of the current split.
    !   splitId  [in] : Integer, split id of the current split.
    !   fa_id    [in] : Integer, node id of the parent node.
    !   p_out [inout] : Real, output probability array allocated previously.
    !
    ! DETAILS
    !   - Initializes the probability matrix for the child nodes.
    !   - Determines the indexes of observations to be used for probability computation
    !     based on the fill type and missingness of features.
    !   - Calls the update_prob subroutine to compute the probabilities for the child nodes.
    ! Added March/2026
    !---------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(in) :: tree
    integer, intent(in) :: nodeId, splitId, fa_id
    real(dp), intent(inout) :: p_out(:, :)

    ! local variables
    integer, allocatable :: idx(:)
    logical, allocatable :: p_zero(:)
    integer :: n_zero, n_prob, feat, j
    type(idx_prob) :: idx_p
    real(dp) :: thr
    real(dp) :: b_lower, b_upper

    ! initialize the probability matrix for the child nodes
    p_out = 0.0_dp

    ! mask for P <= eps
    p_zero = (tree%p(:, fa_id) <= eps)
    n_zero = count(p_zero)

    ! where P_father <= eps
    !   P_left_child = P_father
    !   P_right_child = 0
    ! idx = indexes where P_father > eps
    if (n_zero > 0) then
      where (p_zero) p_out(:, 1) = tree%p(:, fa_id)
      idx = pack([(j, j=1, tree%n_train)], mask=.not. p_zero)
    else
      idx = [(j, j=1, tree%n_train)]
    end if
    n_prob = size(idx)

    feat = tree%nodes(nodeid)%best(splitid)%feat
    thr = tree%nodes(nodeid)%best(splitid)%thr

    ! set the indexes in idx_p according to the fill type
    if (tree%any_na) then
      call get_idx_prob(idx_p, tree%na, tree%nodes(nodeid)%best(splitid)%n_imp, &
        tree%nodes(nodeid)%best(splitid)%imp_feat, n_prob, idx, tree%fill, feat)
    else
      ! If there are no NA's, no need to save other indexes
      idx_p%n_miss_f = 0
      idx_p%n_miss_any = 0
      idx_p%idx = idx
    end if

    call get_feature_bounds(nodeId, feat, tree%nodes_info, tree%thresholds, b_lower, b_upper)

    ! update the probability for the child node
    call update_prob(tree%dist, idx_p, tree%x_train(:, feat), thr, &
                    b_lower, b_upper, tree%sigma(feat), tree%fill, tree%p(:, fa_id), p_out)
  end subroutine update_pmatrix

  subroutine search_from_center(tree, nodeid, try, best, helper, b_lower, b_upper, idx_p, prob)
    !---------------------------------------------------------------------------------
    ! Perform the search for the best split from the middle to the end of the vector.
    ! The direction of the search depends on the 'direc' argument.
    !
    ! ARGUMENTS
    !   tree    [in]    :: tree_model, tree model.
    !   nodeid  [in]    :: Integer, node id of the current split.
    !   try     [inout] :: info_split, information for the split candidate
    !   best    [inout] :: info_split, information on the best splits
    !   helper  [inout] :: h_data, helper variables.
    !   b_lower [in]    :: Real, lower bound for current feature
    !   b_upper [in]    :: Real, upper bound for current feature
    !   idx_p   [in]    :: idx_prob, information on missing values
    !   prob    [in]    :: Real, vector of probabilites for the father node
    !
    ! DETAILS
    !  - Loop over all threholds and compute sum_left and sum_right using Xf and y_cs
    !  - assign_missing: assign missing cases to either node based on a secondary proxy and
    !    update sum_left and sum_right.
    !  - compute the proxy_mse score and compare update the best list for the node.
    !  - Update the probability matrix and compute Prop > min_prob. Ignore threshols that
    !    lead to probabilities too low.
    !
    ! Added July/2025
    ! Last update: May/2026
    !  - replaced several arguments with the tree object to simplify the argument list.
    !--------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: nodeid
    type(info_split), intent(inout) :: try
    type(info_split), intent(inout) :: best(:)
    type(h_data), intent(inout) :: helper
    type(idx_prob), intent(in) :: idx_p
    real(dp), intent(in) :: b_lower, b_upper
    real(dp), intent(in) :: prob(:)
    integer :: i, ii
    logical :: fail
    real(dp) :: y_sum(2), perc_comp(2), multi
    real(dp), allocatable :: p(:, :)

    multi = 1.0_dp / real(tree%n_train, dp)
    allocate(p(tree%n_train, 2))
    i = helper%n_start
    ! stop if gets to any end of the vector
    loop: do
      ! increment i based on the search direction (left or right)
      i = i + helper%direction
      if (i < 1 .or. i > tree%nodes(nodeid)%thresholds(try%feat)%nt) return

      ! save the current threshold
      try%thr = tree%nodes(nodeid)%thresholds(try%feat)%thr(i)

      ! Find the index ii such that x(ii) > threshold
      ! ii > i, always. The equality ii = i + 1 only holds if x contains only
      ! unique values that are not close to each other.
      ii = i
      do while (ii < helper%n_comp)
        if (helper%x(ii + 1) >= try%thr) exit
        ii = ii + 1
      end do

      ! STEP1: Compute left sum and right sum to compute the proxy score using only the observed values
      y_sum(1) = helper%y_cs(ii)
      y_sum(2) = helper%y_cs(helper%n_comp) - y_sum(1)
      try%n_l = ii
      try%n_r = helper%n_comp - try%n_l

      ! Check if the node can be splitted using min_prob (complete analysis)
      !  - compute the probability for the left and right nodes
      !  - compute the proportion of probabilities higher than a threshold
      p = 0.0_dp
      call update_prob(tree%dist, idx_p, tree%x_train(:, try%feat), try%thr, &
                       b_lower, b_upper, tree%sigma(try%feat), tree%fill, prob, p)
      perc_comp = count(p > tree%min_prob, dim=1) * multi

      ! When the first split fails (min_obs or perc_comp), any split from here will also fail
      if (.not. minval(perc_comp) > tree%min_prop) return

      ! If a split can be made, save the regions for the new nodes
      try%regid(helper%idx_c(1:ii)) = left_id                  ! left node
      try%regid(helper%idx_c(ii + 1:helper%n_comp)) = right_id ! right node

      ! STEP2: update left and right sum using missing indexes
      ! Assign Xmiss to some reg using a proxy
      if (helper%n_miss_f > 0) then
        call assign_missing(try, helper%n_miss_f, helper%y_m, y_sum, helper%idx_m, tree%crit)
      end if

      ! compute the proxy score and update the best list
      try%score = y_sum(1)**2 / try%n_l + y_sum(2)**2 / try%n_r

      if (tree%nodes(nodeid)%n_cand_found < tree%n_cand) then
        best(tree%nodes(nodeid)%n_cand_found + 1) = try
       tree%nodes(nodeid)%n_cand_found = tree%nodes(nodeid)%n_cand_found + 1
        cycle
      end if

      call update_best_list(tree%n_cand, best, try, fail)
      if (.not. fail) tree%nodes(nodeid)%n_cand_found = tree%nodes(nodeid)%n_cand_found + 1
    end do loop
  end subroutine search_from_center

  subroutine find_node_splits(tree, nodeid, prob)
    !---------------------------------------------------------------------------------------------
    ! Finds the best candidate splits for a node based on a proxy improvement measure.
    !
    ! ARGUMENTS
    !   tree   [inout] : tree_model, tree model.
    !   nodeid [in]    : Integer, node id of the current split.
    !   prob   [in]    : Real, vector of probabilites for the candidate node
    !
    ! DETAILS
    !  - Initializes P for the two new nodes and identify the indexes where P(father) <= eps.
    !  - Loop over the features to find the best split.
    !    Fore each feature:
    !     - Identify complete cases, sort Xf and Y, compute cumsum(Y) to simplify the search.
    !     - search_from_the_center: Loop over the thresholds (from the center to the borders)
    !       to find the best n_cand candidates
    !
    ! Added July/2025
    ! Last update: May/2026
    !  - removed state and added nodeid
    !---------------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: nodeid
    real(dp), intent(in) :: prob(:)

    ! Local Variables
    integer :: n_zero
    logical, allocatable :: p_zero(:)
    integer :: j, f, n_obs_node, n_prob
    integer :: c_c, c_m
    integer, allocatable :: idx_node(:)
    integer, allocatable :: n_miss(:), idx(:)
    type(info_split), allocatable :: best(:)
    type(info_split) :: try
    type(h_data) :: helper
    type(idx_prob) :: idx_p
    real(dp) :: b_lower, b_upper

    ! Number of observations and their indices in the current node
    ! Each child node will have more than min_obs observations (previously computed)
    n_obs_node = tree%nodes(nodeid)%n_obs_node
    idx_node = tree%nodes(nodeid)%idx  ! global indexes
    allocate (n_miss(tree%n_feat))
    n_miss = count(tree%na(idx_node, :), dim=1) ! number of missing cases in the node

    ! Initialize candidate search for the current node
    tree%nodes(nodeid)%n_cand_found = 0
    allocate(best(tree%n_cand))
    best%score = 0.0_dp
    allocate(try%regid(n_obs_node), source = 0)
    allocate(try%imp_feat(tree%n_feat))

    ! To preserve probability, if the father's node is near zero,
    ! all the probability is kept in the left node.
    allocate(p_zero(tree%n_train))
    p_zero = prob <= eps
    n_zero = count(p_zero)
    if (n_zero > 0) then
      ! Global: pack non zero indexes (need update to identify NA entries)
      idx = pack([(j, j=1, tree%n_train)], mask=.not. p_zero)
    else
      ! Global: Non zero indexes (need update to identify NA entries)
      idx = [(j, j=1, tree%n_train)]
    end if
    n_prob = size(idx)

    ! Pre-allocate matrices to their maximum capacity only ONCE
    ! This avoids expensive OS heap allocation calls during the feature loop
    allocate(helper%idx_c(n_obs_node))
    allocate(helper%x(n_obs_node))
    allocate(helper%y_cs(n_obs_node))
    allocate(helper%idx_m(n_obs_node))
    allocate(helper%y_m(n_obs_node))

    ! Loop over all features to find best splits
    do f = 1, tree%n_feat

      ! Skip feature if it has no valid thresholds (pre-calculated)
      ! if nt > 0, then n_obs > min_obs for each child node (pre-calculated)
      if (tree%nodes(nodeid)%thresholds(f)%nt == 0) cycle

      ! Pre-process complete cases:
      helper%n_miss_f = n_miss(f)
      helper%n_comp = n_obs_node - helper%n_miss_f

      ! Unified single-pass loop replaces three packs, array constructors and slicing evaluations
      c_c = 1
      c_m = 1
      do j = 1, n_obs_node
        if (.not. tree%na(idx_node(j), f)) then
          helper%idx_c(c_c) = j ! Save local index
          helper%x(c_c) = tree%x_train(idx_node(j), f)
          helper%y_cs(c_c) = tree%y_train(idx_node(j))
          c_c = c_c + 1
        else
          helper%idx_m(c_m) = j ! Save local index
          helper%y_m(c_m) = tree%y_train(idx_node(j))
          c_m = c_m + 1
        end if
      end do

      call sort_xy(helper%n_comp, helper%x, helper%y_cs, helper%idx_c)
      do j = 2, helper%n_comp
        helper%y_cs(j) = helper%y_cs(j - 1) + helper%y_cs(j)
      end do

      ! Set current feature for the candidate split
      try%feat = f
      try%n_imp = tree%nodes(nodeid)%n_imp
      try%imp_feat = tree%nodes(nodeid)%imp_feat

      ! Select only the important features to update P
      ! Only needs update if f was not an important feature
      if (try%n_imp == 0 .or. all(try%imp_feat(1:try%n_imp) /= f)) then
        call update_imp_features(try%feat, try%n_imp, try%imp_feat)
      end if

      ! Initialize idx_p to a null state to prevent uninitialized use

      ! set the indexes in idx_p according to the fill type
      if (tree%any_na) then
        call get_idx_prob(idx_p, tree%na, try%n_imp, try%imp_feat, n_prob, idx, tree%fill, f)
      else
        ! If there are no NA's, no need to save other indexes
        idx_p%n_miss_f = 0
        idx_p%n_miss_any = 0
        idx_p%idx = idx
      end if

      call get_feature_bounds(nodeid, f, tree%nodes_info, tree%thresholds, b_lower, b_upper)

      ! Loop in thresholds: starts from the center and goes to the ends (left first)
      ! Stop going in the current direction if a split leads to a region with probability too small

      if (tree%nodes(nodeid)%thresholds(f)%nt > 1) then  ! More likely
        ! go left (nt/2 to 1)
        helper%n_start = tree%nodes(nodeid)%thresholds(f)%nt / 2 + 1
        helper%direction = -1
        call search_from_center(tree, nodeid, try, best, helper, b_lower, b_upper, idx_p, prob)
        ! go right (nt/2 + 1 to nt)
        helper%n_start = tree%nodes(nodeid)%thresholds(f)%nt / 2
        helper%direction = 1
        call search_from_center(tree, nodeid, try, best, helper, b_lower, b_upper, idx_p, prob)
      else
        ! go right (only 1 to test)
        helper%n_start = 0
        helper%direction = 1
        call search_from_center(tree, nodeid, try, best, helper, b_lower, b_upper, idx_p, prob)
      end if
    end do  ! end of loop in features

    ! check if any split can be attempted.
    ! - if n_cand_found == 0, then no split is possible and the node becomes terminal
    if (tree%nodes(nodeid)%n_cand_found == 0) return

    ! set the best candidates and update the node updating status
    tree%nodes(nodeid)%n_cand_found = min(tree%nodes(nodeid)%n_cand_found, tree%n_cand)
    tree%nodes(nodeid)%best = best(1:tree%nodes(nodeid)%n_cand_found)
    tree%nodes(nodeid)%update = .false.
  end subroutine find_node_splits

  subroutine update_mse(tree, nodeid, splitid, gammahat, yhat, mse, p_work)
    !-------------------------------------------------------------------------------------------------
    ! Efficiently updates the coefficients gammahat, the predictions yhat, and the proxy mean
    ! squared error (proxy_mse) for a given split in the tree structure.
    !
    ! ARGUMENTS
    !   tree     [in]    : tree_model, current tree structure.
    !   nodeid   [in]    : Integer, index of the node to split.
    !   splitid  [in]    : Integer, index of the split to apply.
    !   gammahat [out]   : Real, the updated coefficients.
    !   yhat     [out]   : Real, the updated predictions.
    !   mse      [out]   : Real, the updated proxy mse.
    !   p_work   [inout] : Real, working array pre-allocated to avoid overhead.
    !
    ! Added: July/2025
    !-------------------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: nodeid, splitid
    real(dp), allocatable, intent(out) :: gammahat(:)
    real(dp), intent(out) :: yhat(:)
    real(dp), intent(out) :: mse
    real(dp), intent(out) :: p_work(:, :)
    integer :: fa_id

    ! Finds which column of P corresponds to the node being split
    if (tree%n_tn > 1) then
      fa_id = findloc(tree%tn_id(1:tree%n_tn) == nodeid, value=.true., dim=1)
    else
      fa_id = 1
    end if

    ! Fill P with old values, except the ones from current region (jth = fa_id)
    !  - Previous P has n_tn columns.
    !    The columns are indexed by terminal nodes in increasing order
    !  - In the new P, the first n_tn - 1 columns correspond to the old terminal nodes
    !    (which remain the same) and the last two columns correspond to the new regions.
    !    After the jth column, the columns will be shifted one position to the left.
    !  - If P has only one column, there are no columns to shift.
    if (tree%n_tn > 1) then
      if (fa_id > 1) p_work(:, 1:fa_id - 1) = tree%p(:, 1:fa_id - 1)
      if (fa_id < tree%n_tn) p_work(:, fa_id:tree%n_tn - 1) = tree%p(:, fa_id + 1:tree%n_tn)
    end if
    call update_pmatrix(tree, nodeid, splitid, fa_id, p_work(:, tree%n_tn:tree%n_tn + 1))

    ! Update gammahat and the proxy mse value
    gammahat = lsquare(tree, p_work)
    yhat = matmul(p_work, gammahat)
    mse = sum((tree%y_train - yhat)**2)
  end subroutine update_mse

  subroutine update_best_tree_net(tree, split_id)
    !-----------------------------------------------------------------------------------
    ! Updates the tree structure with the best split information.
    !
    ! ARGUMENTS
    !   tree      [inout] : tree_model, current tree structure.
    !   split_id  [in]    : Integer(2), IDs of the nodes to split.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: split_id(2)
    integer :: nodeid, splitid
    integer :: newnode(2), n_nodes, n_tn, feat, fa_id
    real(dp) :: thr
    type(info_split) :: best(1)
    real(dp), allocatable :: new_p(:, :)

    nodeid = split_id(1)
    splitid = split_id(2)
    n_nodes = tree%n_nodes
    n_tn = tree%n_tn
    feat = tree%nodes(nodeid)%best(splitid)%feat
    thr = tree%nodes(nodeid)%best(splitid)%thr

    ! Update P, P_zero, n_zero, gamma, yhat and tn_id before expanding the net
    ! Find the position of the father node in tn_id
    if (n_tn > 1) then
      fa_id = findloc(tree%tn_id(1:n_tn) == nodeid, value=.true., dim=1)
    else
      fa_id = 1
    end if

    ! 1. Calculate new probabilities safely into new_p to avoid aliasing with tree%p
    allocate(new_p(tree%n_train, 2))
    call update_pmatrix(tree, nodeid, splitid, fa_id, new_p)

    if (n_tn > 1) then
      if (fa_id < n_tn) then
         tree%tn_id(fa_id:n_tn - 1) = tree%tn_id(fa_id + 1:n_tn)
         tree%p(:, fa_id:n_tn - 1) = tree%p(:, fa_id + 1:n_tn)
      end if
      tree%tn_id(n_tn:n_tn + 1) = [n_nodes + 1, n_nodes + 2]
      tree%p(:, n_tn:n_tn + 1) = new_p
    else
      tree%tn_id(1:2) = [2, 3]
      tree%p(:, 1:2) = new_p
    end if

    tree%n_nodes = n_nodes + 2
    tree%n_tn = n_tn + 1

    ! Update the old node information
    tree%nodes_info(nodeid, idx_ni_terminal) = 0
    tree%nodes_info(nodeid, idx_ni_feature) = feat
    tree%thresholds(nodeid) = thr
    tree%nodes(nodeid)%split = .false.

    ! Update the new nodes basic information
    newnode = [1, 2] + n_nodes
    tree%nodes_info(newnode, idx_ni_id) = [left_id + n_nodes, right_id + n_nodes]
    tree%nodes_info(newnode, idx_ni_terminal) = 1
    tree%nodes_info(newnode, idx_ni_father) = nodeid
    tree%nodes_info(newnode, idx_ni_depth) = tree%nodes_info(nodeid, idx_ni_depth) + 1
    tree%nodes_info(newnode, idx_ni_feature) = feat
    tree%thresholds(newnode) = thr
    tree%nodes(newnode)%split = .true.

    ! save only the best candidate information for the parent node
    tree%nodes(nodeid)%n_cand_found = 1
    best(1) = tree%nodes(nodeid)%best(splitid)
    tree%nodes(nodeid)%best = best

    ! update the region variable in the net
    tree%region(tree%nodes(nodeid)%idx) = best(1)%regid + n_nodes

    ! Initializing the state variable. New nodes always need update
    ! (other variables will be set after checking basic stopping criterias, if needed)
    tree%nodes(newnode(1))%update = .true.
    tree%nodes(newnode(2))%update = .true.

    ! reset counter
    tree%n_cand_found = 0
  end subroutine update_best_tree_net

  subroutine split_full_analysis(tree, best_id, n_cand_found)
    !-----------------------------------------------------------------------------
    ! Performs a full analysis of the best candidate splits.
    ! This updates the tree structure and selects the best split based on the
    ! lowest mean squared error (mse).
    ! Uses the cp criterion to check if the split is worth it.
    !
    ! ARGUMENTS
    !   tree         [inout] : tree_model, tree object to update.
    !   best_id      [in]    : Integer, IDs of the best candidates.
    !   n_cand_found [in]    : Integer, number of candidates to check.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: n_cand_found
    integer, intent(in) :: best_id(:, :)
    real(dp) :: mse_temp, mse_best
    real(dp), allocatable :: gamma_temp(:), gamma_best(:)
    real(dp), allocatable :: yhat_temp(:), yhat_best(:)
    real(dp), allocatable :: p_work(:, :)
    integer :: k, nodeid, splitid, split_id(2)

    allocate(yhat_temp(tree%n_train), source = 0.0_dp)
    allocate(yhat_best(tree%n_train), source = 0.0_dp)
    allocate(p_work(tree%n_train, tree%n_tn + 1), source = 0.0_dp)
    mse_best = tree%mse_train ! sum of squares of residuals
    split_id = -1

    do k = 1, n_cand_found
      nodeid = best_id(k, 1)
      splitid = best_id(k, 2)

      ! Compute the mse for current candidate and check if the new split improves the mse
      ! Recompute gammahat and yhat for the best tree at the end
      call update_mse(tree, nodeid, splitid, gamma_temp, yhat_temp, mse_temp, p_work)

      if (mse_temp < mse_best) then
        ! If it does, update the best net and mse
        ! and save the best split information
        mse_best = mse_temp
        gamma_best = gamma_temp
        yhat_best = yhat_temp
        split_id(1) = nodeid
        split_id(2) = splitid
      end if
    end do

    ! If no better split was found, disables further splitting and returns
    ! No need to update nodes to terminal. The algorithm will stop as there are
    ! no more candidates to split.
    if (split_id(1) < 0) then
      tree%nodes(1:tree%n_nodes)%split = .false.
      return
    end if

    ! If a split was found, check if the split is worth it using cp criterion.
    ! If the split is not worth it, disables further splitting and returns.
    if (mse_best > (1 - tree%cp) * tree%mse_train) then
      tree%nodes(1:tree%n_nodes)%split = .false.
      return
    end if

    ! full update the tree using the best candidate
    tree%gammahat = gamma_best
    tree%yhat_train = yhat_best
    tree%mse_train = mse_best
    call update_best_tree_net(tree, split_id)
  end subroutine split_full_analysis

  subroutine build_tree(tree)
    !-------------------------------------------------------------------------------------
    !  Builds the probabilistic tree using the provided data and parameters.
    !  WARNING: Before calling build_tree, the tree object must be initialized using
    !    - set_tree_model: to initialize the tree_data/tree_ctrl/etc variables
    !    - start_root: to initialize the root node
    !
    ! ARGUMENTS
    !   tree [inout] : tree_model, tree object to build.
    !
    ! DETAILS
    !  Stopping criteria used
    !   - number of terminal nodes (max_tn)
    !   - number of observations in the final nodes (min_obs)
    !   - node depth (max_d)
    !   - percentage of probabilities (min_prop) higher than a threshold (min_prob)
    !   - reduction in mse (cp)
    !
    !  The top k method is used to choose the node to split and which variable and
    !  threshold to use. The basic steps are
    !   - STEP 1: select the best n_cand candidates to split (by node or globally)
    !     based on a proxy mse criterion.
    !   - STEP 2: performs the full analysis to find the best split among the candidates.
    !   - STEP 3: updates the tree using the best candidate (if any).
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, allocatable :: do_split(:)
    integer :: n_cand_found
    integer :: i, n_split, nodeid, n_best, col_id
    integer, allocatable :: best_id(:, :)
    real(dp), allocatable :: best_score(:)
    logical :: first, fail
    integer :: max_best

    ! Pré-aloca a capacidade máxima possível fora do laço para evitar alocações repetidas
    max_best = tree%n_cand
    if (tree%by_node) max_best = max_best * tree%max_tn
    allocate(best_id(max_best, 2), source = 0)
    allocate(best_score(max_best), source = 0.0_dp)

    ! Main loop to create divisions.
    !  - Uses a stopping criteria based on the number of terminal nodes (max_tn)
    grow_loop: do while ((tree%n_tn < tree%max_tn))

      ! No need to check when n_tn = 1 (always ok)
      ! If n_tn > 1 update do_split and n_split using basic stopping criteria
      !  - criteria used in this step: depth, min_prop, min_prob and max_d
      !    (only needs to check the last two childs)
      if (tree%n_tn == 1) then
        n_split = 1
        do_split = [1]
      else
        ! Updates the list of splittable nodes
        call update_split_candidates(tree, n_split, do_split)
        ! If no further splits are possible we are done growing the three
        if (n_split == 0) then
          if (tree%printinfo >= log_detailed) then
            ! Debug: Print node splitting info
            call labelpr("    No more candidates to split", -1)
            call labelpr(" ", -1)
          end if
          return
        end if

        ! If n_tn = 1 thresholds and NA masks were already computed during initialization
        ! Otherwhise, update the thresholds information and NA masks for the nodes created in the last iteration
        ! (to avoid unnecessary updates child nodes were not updated when the parent node was splitted)
        call update_node_state(tree, fail)    ! left and right child

        ! if the last two nodes are the only splittable nodes and there are no more
        ! thresholds then the search for a new tree is over (no need to update to terminal)
        if (fail) then
           if (n_split <= 2 .and. do_split(1) >= tree%n_nodes - 1) return
        end if
      end if

      ! Atempting a new split.
      !  - If n_tn = 1, finds the best candidates
      !  - if n_tn > 1 loop over the new nodes (if splittable) and update the list of candidates
      if (tree%printinfo >= log_detailed) then
        ! Debug: Print node splitting info
        call labelpr("**************************************************", -1)
        call labelpr("    Attempting a new split ", -1)
        call labelpr("**************************************************", -1)
        call labelpr("Current status", -1)
        call intpr1("    Nodes (n_nodes):", -1, tree%n_nodes)
        call intpr1("    Terminal nodes (n_tn):", -1, tree%n_tn)
        call intpr("    Candidates to split:", -1, do_split, n_split)
        call labelpr(" ", -1)
      end if

      ! Setting variables used int the search for the best candidates
      ! - by_node: controls how the search is done,
      !     if .true., n_cand candidates for each node
      !     if .false., n_cand candidates globally
      ! - first: controls if best_id is to be initialized or updated
      ! - best_id and best_score: informations about the best candidates

      n_best = tree%n_cand          ! at least 1 candidate is needed
      if (tree%by_node) n_best = n_best * n_split
      first = .true.

      ! step 1: uses proxy mse score
      !   - Loop over all nodes to find the n_cand best candidates (by node or globally).
      !   - Eliminate the nodes that do not have enough observations (min_obs criterion).
      !     Missing values are not considered when applying the min_obs criterion.
      !   - net starts with n_cand_found = 0.
      !   - old nodes already have n_cand_found >= 0.
      n_cand_found = 0
      loop_search: do i = 1, n_split
        ! the current candidate for splitting
        nodeid = do_split(i)
        if (.not. tree%nodes(nodeid)%split) cycle

        ! Check if the the node needs update
        ! If update = .false. then the best candidates for this node remain the same
        if (.not. tree%nodes(nodeid)%update) then
          ! Update the best candidates list
          call update_best_list_global(n_best, best_id(1:n_best, :), best_score(1:n_best), &
            nodeid, tree%nodes(nodeid)%n_cand_found, tree%nodes(nodeid)%best, tree%by_node, &
            first, n_cand_found)
          cycle
        end if

        if (tree%printinfo >= log_detailed) then
          ! Debug: Print node splitting info
          call intpr1("    Searching candidates for node Id:", -1, nodeid)
          call intpr1("    Node observations (n_obs_node):", -1, tree%nodes(nodeid)%n_obs_node)
          if (tree%printinfo >= log_debug_deep) then
            call intpr("    Indexes", -1, tree%nodes(nodeid)%idx, tree%nodes(nodeid)%n_obs_node)
          end if
          call labelpr(" ", -1)
        end if

        ! loop over all variables to find the best n_cand for the current node
        col_id = findloc(tree%tn_id(1:tree%n_tn) == nodeid, value=.true., dim=1)
      call find_node_splits(tree, nodeid, tree%p(:, col_id))

        ! disabling further split attempts for nodes that do not have enough data.
        if (tree%nodes(nodeid)%n_cand_found == 0) then
          call update_node_to_terminal(tree%nodes(nodeid))
          cycle
        end if

        ! if at least one candidate was found, update the best candidates list
        call update_best_list_global( &
          n_best, best_id(1:n_best, :), best_score(1:n_best), nodeid, tree%nodes(nodeid)%n_cand_found, &
          tree%nodes(nodeid)%best, tree%by_node, first, n_cand_found)
      end do loop_search

      ! if no candidates were found, there is no split possible
      ! and all splittable nodes were already set as terminal
      if (n_cand_found == 0) return

      ! step 2: perform the full analysis to find the best split
      !   - for each candidate, update the tree nodes with the best split
      !   - choose the best split among the candidates
      !   - uses the cp to check if the split is worth it
      !   - if the split is not worth it, all nodes are set as terminal
      ! On exit, returns the updated tree if a good candidate was found
      n_cand_found = min(n_cand_found, n_best)
      call split_full_analysis(tree, best_id(1:n_cand_found, :), n_cand_found)
    end do grow_loop
  end subroutine build_tree

  !===============================================================================
  ! STEP 4: Predict using validation/test data
  !===============================================================================
  function prob_rgivenx_no_na(argsd, n_test, x, b_lower, b_upper, sigma) result(prob)
    !-------------------------------------------------------------------------------------
    ! Compute the probability P(R | X) = Psi(X, R, sigma) when there are no missing values
    ! and only one feature. Uses the vectorized form of pdist
    !
    ! ARGUMENTS
    !   argsd   [in] : argsDist, distribution related parameters.
    !   n_test  [in] : Integer, number of test observations.
    !   X       [in] : Real, vector of feature values (n_test).
    !   b_lower [in] : Real, lower region bound.
    !   b_upper [in] : Real, upper region bound.
    !   sigma   [in] : Real, standard deviations for each feature.
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n_test
    real(dp), intent(in) :: x(n_test)
    real(dp), intent(in) :: b_lower, b_upper
    real(dp), intent(in) :: sigma
    real(dp), allocatable :: prob(:)

    allocate(prob(n_test))
    ! The conditional checks on bounds are structured to handle the most common case (finite bounds) first.
    if (b_lower > neg_inf .and. b_upper < pos_inf) then
      ! Case 1: R = [lower, upper]. (Most likely)
      prob = pdist(argsd, n_test, b_upper, x, sigma) - pdist(argsd, n_test, b_lower, x, sigma)
    else if (b_lower <= neg_inf) then
      ! Case 2: Lower bound is -inf
      if (b_upper < pos_inf) then
        ! Subcase 2.1: R = (-inf, upper].
        prob = pdist(argsd, n_test, b_upper, x, sigma)
      else
        ! Subcase 2.2: R = (-inf, +inf). (Very unlikely)
        prob = 1.0_dp
      end if
    else
      ! Case 3: R = [lower, +inf).
      prob = 1.0_dp - pdist(argsd, n_test, b_lower, x, sigma)
    end if
  end function prob_rgivenx_no_na

  subroutine prob_rgivenx_sf(n_test, x, na, n_tn, depth, bounds, p, sigma, argsd)
    implicit none
    integer, intent(in) :: n_test, n_tn
    integer, intent(in) :: depth(n_tn)
    real(dp), intent(in) :: x(n_test), bounds(n_tn, 2), sigma
    logical, intent(in) :: na(n_test)
    real(dp), intent(inout) :: p(n_test, n_tn)
    type(argsdist), intent(in) :: argsd
    integer :: i, k, n_comp
    logical :: is_complete
    integer, allocatable :: idx(:)

    ! check if the data has missing values
    is_complete = .not. any(na)

    ! DIRECT APPROACH: Optimization for complete data and single feature
    if (is_complete) then
      do k = 1, n_tn
        p(:, k) = prob_rgivenx_no_na(argsd, n_test, x, bounds(k, 1), bounds(k, 2), sigma)
      end do
      return
    end if

    ! (a) Process all entries that do not have missing values.
    idx = pack([(i, i=1, n_test)], mask=.not. na)
    n_comp = size(idx)
    if (n_comp > 0) then
      ! DIRECT APPROACH: Optimization for complete data and single feature.
      ! The maximum number of calls to pdist is 2 * n_tn * n_test
      do k = 1, n_tn
        p(idx, k) = prob_rgivenx_no_na(argsd, n_comp, x(idx), bounds(k, 1), bounds(k, 2), sigma)
      end do
    end if

    ! (b) If no entries with missing values remain, exit.
    if (n_comp == 0) then
      idx = [(i, i=1, n_test)]
    else
      idx = pack([(i, i=1, n_test)], mask=na)
    end if

    ! (c) Process missing entries hierarchically. Since there exist only one feature,
    !     P(Left) = P(Right) = P(Father)/2 = 1/2^depth
    !     This is calculated directly for each terminal node.
    do k = 1, n_tn
      if (depth(k) < max_prob_depth) then
        p(idx, k) = scale(1.0_dp, -depth(k))
      else
        p(idx, k) = 0.0_dp
      end if
    end do
  end subroutine prob_rgivenx_sf

  function prob_rgivenx(argsd, n_test, n_feat, x, sigma, fill, n_tn, tn, nodes_info, thresholds) result(p)
    !-----------------------------------------------------------------------------
    ! Computes probability matrix P for new data with missing values.
    ! The calculation is done hierarchically to mimic the procedure followed
    ! during the tree building process.
    !
    ! ARGUMENTS
    !   argsd      [in] : argsDist, distribution related parameters.
    !   n_test     [in] : Integer, number of test observations.
    !   n_feat     [in] : Integer, number of features.
    !   X          [in] : Real, data matrix (n_test x n_feat).
    !   sigma      [in] : Real, vector of sigmas (n_feat).
    !   fill       [in] : Integer, fill strategy for missing values.
    !   n_tn       [in] : Integer, number of terminal nodes.
    !   tn         [in] : Integer, terminal nodes index.
    !   nodes_info [in] : Integer, information on nodes. (2*n_tn - 1, :)
    !                     (id* optional, isTerminal, father, depth, feature)
    !   thresholds [in] : Real, thresholds for the nodes.
    !   p          [out]: Real, the resulting probability matrix
    !                     of size (n_test x n_tn).
    !
    ! DETAILS
    !   This function uses two main strategies based on the number of features:
    !   1. n_feat = 1: A hybrid approach is used.
    !      - Direct calculation for observations without missing values.
    !      - Handles missing values using P(child) = 1/2^depth
    !   2. n_feat > 1: The hierarchical approach is used for all observations.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n_test, n_feat, n_tn, fill
    integer, intent(in) :: tn(:)
    integer, intent(in) :: nodes_info(:, :)
    real(dp), intent(in) :: x(:, :), sigma(:)
    real(dp), intent(in) :: thresholds(:)
    real(dp), allocatable :: p(:, :)

    ! Local variables
    integer :: i, k, n_nodes, fa_id, f, offset, ncol_ni
    integer :: n_imp, n_prob, n_zero
    integer :: imp_feat(n_feat)
    integer, allocatable :: idx(:), n_miss(:)
    logical :: is_complete
    logical, allocatable :: na(:, :), p_zero(:)
    real(dp), allocatable :: p_matrix(:, :)
    real(dp), allocatable :: bd1(:, :)
    real(dp) :: b_lower, b_upper
    type(idx_prob) :: idx_p

    allocate(p(n_test, n_tn))
    allocate(na(n_test, n_feat))
    allocate(p_zero(n_test))
    allocate(p_matrix(n_test, 2 * n_tn - 1))
    allocate(n_miss(n_feat))
    allocate(bd1(n_tn, 2))

    ! Offset mapping based on ncol_ni: 0 if 5 columns, 1 if 4 columns
    ncol_ni = size(nodes_info, 2)
    offset = 5 - ncol_ni

    ! Early exit: if only one terminal node, probability is 1 for all observations.
    if (n_tn == 1) then
      p(:, 1) = 1.0_dp
      return
    end if

    ! Total number of nodes in the tree
    n_nodes = 2 * n_tn - 1

    ! Identify missing values in the data matrix
    na = isnan(x)
    p_matrix = 0.0_dp

    if (n_feat > 1) then
      ! Number of missing values by feature
      n_miss = count(na, dim=1)
      ! Features used to build the tree
      n_imp = count([(any(i == nodes_info(:, idx_ni_feature - offset)), i=1, n_feat)])
      if (n_imp > 0) then
        imp_feat(1:n_imp) = pack([(i, i=1, n_feat)], &
          mask=[(any(i == nodes_info(:, idx_ni_feature - offset)), i=1, n_feat)])
      else
        n_imp = 1
        imp_feat(1) = 1
      end if
      ! mask for complete dataset
      is_complete = all(n_miss(imp_feat(1:n_imp)) == 0)
    else
      n_imp = 1
      imp_feat(1) = 1
    end if

    ! Case 1: Single important feature.
    !  - find the indexes for complete data and use the direct approach.
    !  - find the indexes for missing data and use the hierarchical approach (if any)
    if (n_imp == 1) then
      f = imp_feat(1)
      do i = 1, n_tn
        call get_feature_bounds(tn(i), f, nodes_info, thresholds, bd1(i, 1), bd1(i, 2))
      end do
      call prob_rgivenx_sf(n_test, x(:, f), na(:, f), n_tn,&
       nodes_info(tn, idx_ni_depth - offset), bd1, p, sigma(f), argsd)
      return
    end if

    ! Case 2: Multiple features. Process all entries hierarchically.
    ! HIERARCHICAL APPROACH:  The maximum number of calls to pdist is 3 * (n_tn - 1) * n_test.
    ! Loop through sibling pairs to compute probabilities.

    p_matrix(:, 1) = 1.0_dp  ! root node
    do k = 2, n_nodes - 1, 2
      ! find the current father node id and splitting feature
      fa_id = nodes_info(k, idx_ni_father - offset)
      f = nodes_info(fa_id, idx_ni_feature - offset)

      ! To preserve probability, if the father's node is near zero,
      ! all the probability is kept in the left node.
      p_zero = (p_matrix(:, fa_id) <= eps)
      n_zero = count(p_zero)
      if (n_zero > 0) then
        where (p_zero)
          p_matrix(:, k) = p_matrix(:, fa_id)
        end where
        ! pack non zero indexes (needs update to remove NA entries)
        idx = pack([(i, i=1, n_test)], mask=.not. p_zero)
      else
        ! pack non zero indexes (needs update to remove NA entries)
        idx = [(i, i=1, n_test)]
      end if
      n_prob = size(idx)

      call get_node_imp_features(k, nodes_info, imp_feat, n_imp)

      ! if there are no missing values, then only checks the p_zero mask
      if (is_complete) then
        idx_p%n_miss_f = 0
        idx_p%n_miss_any = 0
        idx_p%idx = idx ! p > 0
      else
        ! na gives the mask for missing values while idx gives the indexes of P > eps
        ! Xf could be missing but P(father) <= eps
        call get_idx_prob(idx_p, na, n_imp, imp_feat, n_prob, idx, fill, f)
      end if

      call get_feature_bounds(fa_id, f, nodes_info, thresholds, b_lower, b_upper)

      ! update the probability matrix
      call update_prob(argsd, idx_p, x(:, f), thresholds(fa_id), b_lower, b_upper, &
        sigma(f), fill, p_matrix(:, fa_id), p_matrix(:, k:k + 1))
    end do
    p = p_matrix(:, tn)
  end function prob_rgivenx

  subroutine predict_tree(tree)
    implicit none
    type(tree_model), intent(inout) :: tree
    integer :: n_feat
    real(dp), allocatable :: yhat(:)

    if (tree%n_test == 0) then
      ! If there is no holdout set, use the training MSE as the selection criterion.
      tree%mse_test = tree%mse_train
      return
    end if

    ! Compute predictions and validation/test MSE for the current candidate tree.
    allocate(yhat(tree%n_test))

    ! Set the information to pass to the probability function
    n_feat = tree%n_feat

    !---------------------------------------------
    ! filling the matrix P with the values
    !    Psi(x, Rj, sigma), 1 <= j <= n_tn
    ! where Rj is the j-th terminal node (reg)
    !---------------------------------------------
    tree%p_test = 0.0_dp
    tree%p_test(:, 1:tree%n_tn) = prob_rgivenx(tree%dist, tree%n_test, n_feat, tree%x_test, &
      tree%sigma, tree%fill, tree%n_tn, tree%tn_id(1:tree%n_tn), tree%nodes_info, tree%thresholds)
    yhat = matmul(tree%p_test(:, 1:tree%n_tn), tree%gammahat(1:tree%n_tn))
    tree%yhat_test = yhat
    tree%mse_test = sum((tree%yhat_test - tree%y_test)**2)
  end subroutine predict_tree

  !===============================================================================
  ! STEP 5: Extract the information from the treee object
  !===============================================================================
  subroutine return_tree(tree, n_tn, nodes_info, thresholds, yhat_train, yhat_test, mse, xregion, final)
    !----------------------------------------------------------------------------
    ! Returns the final tree structure to R.
    ! This subroutine copies the values from the final tree structure
    ! to the output vectors and matrices.
    !
    ! ARGUMENTS
    !   tree       [in]    : tree_model, fitted tree object
    !   n_tn       [out]   : Integer, number of terminal nodes in the final tree.
    !   nodes_info [out]   : Integer, matrix with information about the nodes
    !                        in the final tree. (2 * max_tn - 1, 5)
    !   thresholds [inout] : Real, vector of thresholds defining the
    !                         regions in the final tree. (2 * max_tn - 1)
    !   yhat_train [inout] : Real, predicted values for Y (n_train).
    !   yhat_test  [inout] : Real, predicted values for Y (n_test).
    !   mse        [inout] : Real, mean square error for the final tree.
    !   Xregion    [inout] : Integer, matrix with the indices of the features
    !                        used in the nodes (n_train).
    !   final      [in]    : Logical, if .true. the code is ready to return
    !                        to R. Otherwise the tree object needs to be restarted
    !----------------------------------------------------------------------------
    ! Last updated: May/2026
    !  - replaced net with the tree object
    implicit none
    type(tree_model), intent(inout) :: tree
    logical, intent(in) :: final
    integer, intent(out), target :: nodes_info(:, :)
    integer, intent(out) :: n_tn
    integer, intent(out) :: xregion(:)
    real(dp), intent(out), target :: thresholds(:)
    real(dp), intent(out) :: yhat_train(:)
    real(dp), intent(out) :: yhat_test(:)
    real(dp), intent(out) :: mse(3)
    integer :: n_nodes

    ! Updating the number of terminal nodes
    ! number of nodes in the final tree
    n_nodes = tree%n_nodes

    ! Updating the number of terminal nodes and P
    n_tn = tree%n_tn

    ! If the tree is not final, the tree object is being used to store the best candidate tree
    if(.not. final) then
      tree%p_temp = tree%p(:, 1:n_tn)
      tree%gamma_temp = tree%gammahat(1:n_tn)
      if(tree%n_test > 0) tree%p_test_temp = tree%p_test(:, 1:n_tn)
    end if

    ! predicted values and mse
    yhat_train = tree%yhat_train
    mse(1) = tree%mse_train / tree%n_train
    if(tree%n_test > 0) then
      yhat_test = tree%yhat_test
      mse(2) = tree%mse_test / tree%n_test
      mse(3) = (tree%mse_train + tree%mse_test) / tree%n_obs
    else
      mse(2) = 0.0_dp
      mse(3) = mse(1)
    end if

    ! updating the nodes_info
    if (.not. associated(tree%nodes_info, nodes_info)) then
    nodes_info = 0
      nodes_info(1:n_nodes, idx_ni_id) = tree%nodes_info(1:n_nodes, idx_ni_id)
      nodes_info(1:n_nodes, idx_ni_terminal) = tree%nodes_info(1:n_nodes, idx_ni_terminal)
      nodes_info(1:n_nodes, idx_ni_father) = tree%nodes_info(1:n_nodes, idx_ni_father)
      nodes_info(1:n_nodes, idx_ni_depth) = tree%nodes_info(1:n_nodes, idx_ni_depth)
      nodes_info(1:n_nodes, idx_ni_feature) = tree%nodes_info(1:n_nodes, idx_ni_feature)
    end if

    ! updating the thresholds
    if (.not. associated(tree%thresholds, thresholds)) then
      thresholds = 0.0_dp
      thresholds(1:n_nodes) = tree%thresholds(1:n_nodes)
    end if

    ! updating Xregion
    xregion = tree%region
  end subroutine return_tree

  subroutine return_root_only_tree(tree, y_test, keep_full, n_tn, yhat_train, yhat_test,  &
    p_test, mse, nodes_info, thresholds, sigma_best, xregion)
    !----------------------------------------------------------------------------
    ! Handles the degenerate case where the root node cannot be split.
    !
    ! In this case, the final tree has only one terminal node (the root).
    ! This routine fills the outputs consistently for both:
    !   - full-output mode
    !   - light/CV mode
    !
    ! ARGUMENTS
    !   tree        [inout] : tree_model, fitted root-only tree.
    !   y_test      [in]    : Real, full response vector.
    !   keep_full   [in]    : Logical, whether full tree output is required.
    !   n_tn        [out]   : Integer, number of terminal nodes.
    !   yhat_train  [inout] : Real, predictions for training set
    !   yhat_test   [inout] : Real, predictions for test set.
    !   p_test      [inout] : Real, probabilities for test set.
    !   mse         [inout] : Real, mean square error for test set.
    !   nodes_info  [inout] : Integer, node metadata.
    !   thresholds  [inout] : Real, split thresholds.
    !   sigma_best  [out]   : Real, selected sigma vector.
    !   Xregion     [out]   : Integer, assigned region for each training observation (n_train).
    !
    ! Added March/2026
    !----------------------------------------------------------------------------
    implicit none

    ! Input
    type(tree_model), intent(inout) :: tree
    real(dp), intent(in) :: y_test(:)
    logical, intent(in) :: keep_full

    ! Output
    integer, intent(out) :: n_tn
    real(dp), intent(inout) :: mse(3)
    real(dp), intent(inout) :: yhat_train(:)
    real(dp), intent(inout) :: yhat_test(:)
    real(dp), intent(inout) :: p_test(:, :)
    integer, intent(inout), target :: nodes_info(:, :)
    real(dp), intent(inout) :: thresholds(:)
    real(dp), intent(out) :: sigma_best(:)
    integer, intent(inout) :: xregion(:)

    ! Local
    integer :: n_test

    if (tree%printinfo >= log_base) then
      ! Debug: Print node splitting info
      call labelpr("Root node cannot be splitted: No thresholds found.", -1)
      call labelpr(" ", -1)
    end if

    n_tn = 1

    ! Root-only tree: no sigma was effectively selected
    sigma_best = 0.0_dp

    ! Number of test observations
    n_test = tree%n_test

    ! updates yhat_train, P_train and Xregion for the training sample
    if (keep_full) then
      ! Serialize the trivial root-only tree to the output buffers
      call return_tree(tree, n_tn, nodes_info, thresholds, &
       yhat_train, yhat_test, mse, xregion, .true.)
       if(n_test > 0) then
         p_test = 0.0_dp
         p_test(:, 1) = 1.0_dp
       end if
    else
      ! In light mode, only the training MSE is needed
      mse(1) = tree%mse_train / tree%n_train
      if(tree%n_test > 0) then
        mse(2) = sum((y_test - tree%gammahat(1))**2) / tree%n_test
        mse(3) = (tree%mse_train + mse(2)) / tree%n_obs
      else
        mse(2) = 0
        mse(3) = mse(1)
      end If
    end if
  end subroutine return_root_only_tree

end module prtree_main
