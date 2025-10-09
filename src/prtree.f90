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
! Part 2: defines the main interface for the PRTree algorithm,
! including training, prediction, and model evaluation functions.
!
! Build upon the original PRTree implementation by Alisson S. Neimaier.
!
!-------------------------------------------------------------------------------
! This version: July/2025
! By Taiane S. Prass - PPGEst/UFRGS
!-------------------------------------------------------------------------------
module prtree_main
  use prtree_types
  use prtree_misc
  implicit none

  interface return_tree
    ! These subroutines extract components from the tree model
    ! either to return to R or to process the test sample
    ! Added July/2025
    module procedure return_tree_train
    module procedure return_tree_test
  end interface return_tree

contains
  !======================================================
  ! STEP 1: start the tree model
  ! Initialization performed:
  !   => tree%train: COMPLETE
  !   => tree%test : NONE (set after, if needed)
  !   => tree%ctrl : COMPLETE
  !   => tree%net  : NONE (set after)
  !   => tree%dist : COMPLETE
  !======================================================
  pure subroutine set_tree_model(tree, n_obs, n_feat, y, x, int_param, dble_param)
    !----------------------------------------------------------------------------
    ! Allocates and initializes the tree structure and input data.
    ! Must be called once before building the tree.
    ! Initialization performed:
    !   => tree%train: COMPLETE: set all variables
    !   => tree%test : NONE (set after, if needed)
    !   => tree%ctrl : COMPLETE: set all variables
    !   => tree%net  : PARTIAL: allocates region and yhat
    !   => tree%dist : COMPLETE: set all variables
    !
    ! ARGUMENTS
    !   tree       [out] : tree_model, tree object to initialize.
    !   n_obs      [in]  : Integer, number of observations.
    !   n_feat     [in]  : Integer, number of features.
    !   y          [in]  : Real(dp), response vector.
    !   X          [in]  : Real(dp), feature matrix.
    !   int_param  [in]  : Integer(:), parameters to build the tree
    !                      [fill, crit, max_tn, max_d, min_obs, n_cand, by_node, dist_id]
    !   dble_param [in]  : Real(dp), parameters to build the tree
    !                      [min_prop, min_prob, cp, par]
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(out) :: tree
    integer, intent(in) :: n_obs, n_feat
    real(dp), intent(in) :: y(n_obs), x(n_obs, n_feat)
    integer, intent(in) :: int_param(8)
    real(dp), intent(in) :: dble_param(4)

    ! Set input data and parameters (all)
    tree%train%n_obs = n_obs
    tree%train%n_feat = n_feat
    tree%train%fill = int_param(1)
    tree%train%y = y
    tree%train%x = x
    tree%train%sigma = vector(n_feat, source=0.0_dp)
    tree%train%na = isnan(x)
    tree%train%n_miss = count(tree%train%na, dim=1)
    tree%train%any_na = sum(tree%train%n_miss) > 0

    ! Set building criteria (all)
    tree%ctrl%crit = int_param(2)
    tree%ctrl%max_tn = int_param(3)
    tree%ctrl%max_d = int_param(4)
    tree%ctrl%min_obs = int_param(5)
    tree%ctrl%n_cand = int_param(6)
    tree%ctrl%by_node = int_param(7) == 1
    tree%ctrl%min_prop = dble_param(1)
    tree%ctrl%min_prob = dble_param(2)
    tree%ctrl%cp = dble_param(3)

    ! Initialize the argsDist object
    tree%dist%dist_id = int_param(8)
    tree%dist%dist_par = dble_param(4)
  end subroutine set_tree_model

  !===================================================================
  ! STEP 2: Initialize the root node
  ! Initialization performed:
  !   => tree%net: COMPLETE  (2.1: start_root)
  !   => tree%net%nodes: COMPLETE (2.2: start_nodes_root)
  !   => tree%net%nodes%state: PARTIAL (2.3: start_node_state_root)
  !     - %imp_feature: will be initialized DURING split (n_imp = 0)
  !     - %best: will be initialized DURING split. (update = .false.)
  !===================================================================
  pure subroutine start_node_state_root(state, n_obs, n_feat)
    !------------------------------------------------------------------------
    ! Initializes the node state and missing value info for the root node.
    ! Some variables will require initialization/update before use
    !   - imp_feature: will be initialized DURING split (n_imp = 0)
    !   - na_info: will be initialized BEFORE split
    !   - thresholds: will be updated BEFORE split
    !   - best: will be initialized DURING split. (update = .false.)
    !
    ! ARGUMENTS
    !   state  [inout] : node_state, state to initialize.
    !   n_obs  [in]    : Integer, number of observations.
    !   n_feat [in]    : Integer, number of features.
    !
    ! Added July/2025
    !------------------------------------------------------------------------
    implicit none
    type(node_state), intent(out) :: state
    integer, intent(in) :: n_obs, n_feat
    type(info_thr) :: null_thr
    integer :: i

    ! Set the initial status
    state%update = .true.    ! need update
    state%n_cand_found = 0   ! no candidates found at this point
    state%n_obs = n_obs      ! all data starts in the root node
    state%idx = [(i, i=1, n_obs)]
    state%thresholds = vector(n_feat, null_thr) ! update after exiting
    state%n_imp = 0 ! to start the update
    state%imp_feat = [0] ! to start the update
  end subroutine start_node_state_root

  pure subroutine start_nodes_root(nodes, n_obs, n_feat)
    !-------------------------------------------------------------
    !  Allocates and initializes the root node and its bounds.
    !
    ! ARGUMENTS
    !   nodes  [inout] : tree_node(:), array of nodes.
    !   n_obs  [in]    : Integer, number of observations.
    !   n_feat [in]    : Integer, number of features.
    !
    ! Added July/2025
    !-------------------------------------------------------------
    implicit none
    type(tree_node), allocatable, intent(out) :: nodes(:)
    integer, intent(in) :: n_obs, n_feat
    type(tree_node) :: node

    ! Allocate the nodes object and the bounds argument
    !  - In the root node, -inf < X_j < inf
    node%id = 1
    node%isterminal = 1
    node%fathernode = 0
    node%depth = 0
    node%feature = 0
    node%threshold = 0.0_dp
    node%split = .true.
    node%bounds = matrix(n_feat, 2, neg_inf)
    node%bounds(:, 2) = pos_inf         ! update the second column

    ! initializes the state variables
    call start_node_state_root(node%state, n_obs, n_feat)
    nodes = [node] ! ensure that the node object is an array
  end subroutine start_nodes_root

  subroutine start_root(tree)
    !--------------------------------------------------------------
    ! Initializes the tree structure and root node.
    ! Initialization performed:
    !   => tree%net: COMPLETE  (2.1: start_root)
    !   => tree%net%nodes: COMPLETE (2.2: start_nodes_root)
    !   => tree%net%nodes%state: PARTIAL (2.3: start_node_state_root)
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
    integer :: f, n, i
    integer, allocatable :: idx(:)

    ! At start number of nodes = number of terminal nodes
    tree%net%n_nodes = 1
    tree%net%n_tn = 1
    tree%net%n_cand_found = 0

    ! In the root node
    ! - Region = 1 and tn_id = 1
    ! - P = 1 and gammahat = mean(y)
    tree%net%region = vector(tree%train%n_obs, 1)
    tree%net%tn_id = [1]
    tree%net%n_zero = [0] ! (all P > eps)
    tree%net%p_zero = matrix(tree%train%n_obs, 1, source=.false.)
    tree%net%p = matrix(tree%train%n_obs, 1, 1.0_dp)
    tree%net%gammahat = [sum(tree%train%y) / tree%train%n_obs]
    tree%net%yhat = vector(tree%train%n_obs, tree%net%gammahat(1))
    tree%net%mse = sum(tree%train%y**2) - tree%train%n_obs * tree%net%gammahat(1)**2

    ! Allocate the root node and its components
    call start_nodes_root(tree%net%nodes, tree%train%n_obs, tree%train%n_feat)

    ! Loop to compute the thresholds for each feature
    !  - checks if the node can be splitted using min_obs
    if (printinfo >= log_debug) then
      ! Debug: thresholds
      call labelpr("------------------------------------------------------------------", -1)
      call labelpr("    Getting the thresholds for the features at the root node", -1)
      call labelpr("------------------------------------------------------------------", -1)
    end if

    loop_find_thr: do f = 1, tree%train%n_feat
      ! Create mask for non-missing values of feature f
      if (tree%train%n_miss(f) == 0) then
        idx = [(i, i=1, tree%train%n_obs)]
      else
        idx = pack([(i, i=1, tree%train%n_obs)], mask=.not. tree%train%na(:, f))
      end if
      n = tree%train%n_obs - tree%train%n_miss(f)

      ! Skip feature if insufficient non-missing values
      if (n < 2 * tree%ctrl%min_obs) then
        tree%net%nodes(1)%state%thresholds(f)%nt = 0
        cycle
      end if

      ! Find thresholds for current feature (considering only non-missing values)
      call find_thresholds(tree%ctrl%min_obs, n, tree%train%x(idx, f), tree%net%nodes(1)%state%thresholds(f))

      if (printinfo >= log_debug) then
        ! Debug: thresholds found
        call intpr1("    Feature ", -1, f)
        if (tree%net%nodes(1)%state%thresholds(f)%nt == 0) then
          call labelpr("    No thresholds found", -1)
        else
          call intpr1("    Number of thresholds:", -1, &
                      tree%net%nodes(1)%state%thresholds(f)%nt)
        end if
        call labelpr(" ", -1)
      end if
    end do loop_find_thr

    ! If there are missing values, set the NA mask for node(1)%state
    ! Here p > eps for all obserations, so just copy the NA mask
    if (tree%train%any_na) then
      ! intitialize the na_info variable
      select case (tree%train%fill)
      case (0, 1)
        ! Both masks required. The mask for any_na will be updated during split
        n = tree%train%n_feat + 1
        tree%net%nodes(1)%state%na_info = matrix(tree%train%n_obs, n, .false.)
        tree%net%nodes(1)%state%na_info(:, 1:n - 1) = tree%train%na
      case (2)
        ! Only na_f is required
        tree%net%nodes(1)%state%na_info = tree%train%na
      end select
    end if
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
  !      If any_na update the na_info masks for child nodes (NA and P > eps)
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
    type(node_state) :: null_state

    ! disable further splits
    node%split = .false.

    ! Information on splitting is no longer required.
    ! Deallocate the components from state variable
    node%state = null_state
  end subroutine update_node_to_terminal

  pure subroutine update_split_candidates(net, ctrl, n_obs, n_split, do_split)
    !---------------------------------------------------------------------------
    ! Checks which nodes are eligible for splitting based on the depth of
    ! the node and probability criteria. Number of observations are checked
    ! in another step
    !
    ! ARGUMENTS
    !   net      [inout] : tree_net, tree structure.
    !   ctrl     [inout] : tree_ctrl, control parameters.
    !   n_obs    [in]    : Integer, number of observations.
    !   n_split  [inout] : Integer, number of splittable nodes (old vs current).
    !   do_split [inout] : Integer, indexes of splittable nodes (old vs current).
    !
    ! Added July/2025
    !---------------------------------------------------------------------------
    implicit none
    type(tree_net), intent(inout) :: net
    type(tree_ctrl), intent(inout) :: ctrl
    integer, allocatable, intent(inout) :: do_split(:)
    integer, intent(in) :: n_obs
    integer, intent(inout) :: n_split
    real(dp) :: perc_comp
    integer :: j, cut(n_split + 2), n_found, tn_id, nodeid

    ! copy the id of the previous splittable nodes
    cut = [do_split, net%n_nodes - 1, net%n_nodes]

    ! check the current split condition of the node
    !  - split condition can change during the search for thresholds
    !    when the number of observations is evaluated.
    !  - the node is updated to terminal when the split condition changes.

    n_found = 0
    do j = 1, n_split + 2
      nodeid = cut(j)

      if (.not. net%nodes(nodeid)%split) cycle ! node updated to terminal during the process

      ! Check stopping criteria (only for last two child nodes):
      !  - depth and the percentage of probabilities higher than a threshold
      !  - use P > 2 * min_prob as a proxy for the next division.
      if (j > n_split) then
        if (net%nodes(nodeid)%depth >= ctrl%max_d) then
          call update_node_to_terminal(net%nodes(nodeid))
          cycle
        end if

        ! P(child) < a, whenever P(father) < a.
        ! If a <= P(father) < 2*a implies P(child1) < a, whenever P(child2) > a
        tn_id = merge(net%n_tn - 1, net%n_tn, j == n_split + 1)
        perc_comp = count(net%p(:, tn_id) > 2.0_dp * ctrl%min_prob) / real(n_obs, dp)
        if (perc_comp <= ctrl%min_prop) then
          call update_node_to_terminal(net%nodes(nodeid))
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

  subroutine print_message(left, f, n)
    !----------------------------------------------------------------------
    ! Helper subroutine to print debug messages about thresholds found
    ! for a given feature in a node.
    !
    ! ARGUMENTS
    !   left [in] : logical, .true. if left node, .false. if right node.
    !   f    [in] : Integer, feature index.
    !   n    [in] : Integer, number of thresholds found.
    !
    ! Added July/2025
    !----------------------------------------------------------------------
    implicit none
    logical, intent(in) :: left
    integer, intent(in) :: f, n
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

  subroutine update_state_thresholds(old_node, min_obs, n_feat, x, na, left_node, right_node, updated)
    !--------------------------------------------------------------------------------------
    ! Updates thresholds for child nodes after a split operation.
    ! Splits the threshold vector for current variable and recompute for others
    !
    ! ARGUMENTS
    !   old_node    [inout] : type(tree_node), last splitted node.
    !                         state%best is expected to contain only the best split.
    !   min_obs     [in]    : Integer, minimum observations per node.
    !   n_feat      [in]    : Integer, number of features.
    !   X           [in]    : Real(dp), feature values (n_obs x n_feat) for splitted node.
    !   na          [in]    : Logical, missing mask (n_obs x n_feat).
    !   split       [in]    : logical, if .false. does not need update.
    !   left_node   [inout] : type(tree_node), left node.
    !   right_node  [inout] : type(tree_node), right node.
    !   updated     [inout] : logical, if .true., new splits can be attempted
    !
    ! NOTES
    !   Thresholds are allocated during split and only updated here
    !
    !   For current variable:
    !    - Checks min_obs requirement for future splits in both directions
    !    - Preserves original threshold ordering
    !
    ! Added July/2025
    !--------------------------------------------------------------------------------------
    implicit none
    type(tree_node), intent(inout) :: old_node
    integer, intent(in) :: min_obs, n_feat
    real(dp), intent(in) :: x(old_node%state%n_obs, n_feat)
    logical, intent(in) :: na(old_node%state%n_obs, n_feat)
    type(tree_node), intent(inout) :: left_node, right_node
    logical, intent(inout) :: updated(2)
    ! Local variables
    integer :: pos, i, f, feat, n_node
    real(dp) :: thr
    logical :: msk_left(old_node%state%n_obs)
    logical :: msk_right(old_node%state%n_obs)
    logical :: mask_miss(old_node%state%n_obs)
    type(info_thr) :: feat_thr
    integer :: n1, n2, n
    integer, allocatable :: idx(:)

    if (printinfo >= log_debug) then
      ! Debug: thresholds found
      call labelpr("---------------------------------------------", -1)
      call labelpr("    Updating thresholds for the new nodes", -1)
      call labelpr("---------------------------------------------", -1)
    end if

    ! splitting feature and threshold
    feat = old_node%feature
    thr = old_node%threshold

    ! mask for left and right nodes
    msk_left = old_node%state%best(1)%regid == left_id
    msk_right = .not. msk_left

    ! counts the number of thresholds found globally (left and right)
    n1 = 0
    n2 = 0

    if (left_node%split) left_node%state%thresholds%nt = 0
    if (right_node%split) right_node%state%thresholds%nt = 0

    do f = 1, n_feat
      ! current feature will be processed at the end of the loop
      if (f == feat) cycle

      ! mask for non-missing values
      mask_miss = .not. na(:, f)

      ! recompute thresholds for left node
      n = 0         ! local counter
      if (left_node%split) then
        idx = pack([(i, i=1, old_node%state%n_obs)], mask=msk_left .and. mask_miss)
        call find_thresholds(min_obs, size(idx), x(idx, f), left_node%state%thresholds(f))
        n = left_node%state%thresholds(f)%nt
        n1 = n1 + n
      end if
      call print_message(.true., f, n)

      ! recompute thresholds for right node
      n = 0
      if (right_node%split) then
        idx = pack([(i, i=1, old_node%state%n_obs)], mask=msk_right .and. mask_miss)
        call find_thresholds(min_obs, size(idx), x(idx, f), right_node%state%thresholds(f))
        n = right_node%state%thresholds(f)%nt
        n2 = n2 + n
      end if
      call print_message(.false., f, n)
    end do

    ! Finds the position of thr in the thresholds vector
    pos = findloc(old_node%state%thresholds(feat)%thr >= thr, value=.true., dim=1)

    ! mask for non missing values
    mask_miss = .not. na(:, feat)

    ! thresholds information for splitting feature
    feat_thr = old_node%state%thresholds(feat)

    ! LEFT NODE: can only split using thresholds < thr
    ! Since these observations are already in left node (msk = .true.),
    ! we only need to ensure the RIGHT child after split has >= min_obs
    n = 0   ! local counter
    if (left_node%split .and. pos > 1) then ! There exist at least one thresholds < thr
      do i = pos - 1, 1, -1                 ! Check from largest < thr downward
        n_node = count(msk_left .and. mask_miss .and. x(:, feat) > feat_thr%thr(i))
        if (n_node >= min_obs) then
          left_node%state%thresholds(feat)%nt = i
          left_node%state%thresholds(feat)%thr = feat_thr%thr(1:i)
          exit
        end if
      end do
      n = left_node%state%thresholds(feat)%nt
      n1 = n1 + n
    end if
    call print_message(.true., feat, n)

    ! RIGHT NODE: can only split using thresholds > thr
    ! Since these observations are already in right node (msk = .false.),
    ! we only need to ensure the LEFT child after split has >= min_obs
    n = 0  ! local counter
    if (right_node%split .and. pos < feat_thr%nt) then ! There exist at least one thresholds > thr
      do i = pos + 1, feat_thr%nt                      ! Start after the exact split threshold
        n_node = count(msk_right .and. mask_miss .and. x(:, feat) <= feat_thr%thr(i))
        if (n_node >= min_obs) then
          right_node%state%thresholds(feat)%nt = feat_thr%nt - i + 1
          right_node%state%thresholds(feat)%thr = feat_thr%thr(i:feat_thr%nt)
          exit
        end if
      end do
      n = right_node%state%thresholds(feat)%nt
      n2 = n2 + n
    end if
    call print_message(.false., feat, n)

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
    !----------------------------------------------------------------------
    implicit none
    type(tree_model), target, intent(inout) :: tree
    logical, intent(inout) :: fail
    type(tree_node), pointer :: p_node, l_node, r_node
    integer :: nodeid, j, n_tn, n
    logical, allocatable :: msk(:)
    logical :: update_zero(2)
    type(info_thr) :: null_thr
    logical :: update(2)

    ! find the node that originated the two last nodes
    nodeid = tree%net%nodes(tree%net%n_nodes)%fathernode
    p_node => tree%net%nodes(nodeid)               ! parent node
    l_node => tree%net%nodes(tree%net%n_nodes - 1) ! left child
    r_node => tree%net%nodes(tree%net%n_nodes)     ! right child

    if (.not. (l_node%split .or. r_node%split)) then
      call update_node_to_terminal(l_node)
      call update_node_to_terminal(r_node)
      call update_node_to_terminal(p_node)
      fail = .true.
      return
    end if

    msk = p_node%state%best(1)%regid == left_id

    if (l_node%split) then
      l_node%state%n_obs = p_node%state%best(1)%n_l                        ! size
      l_node%state%idx = pack(p_node%state%idx, mask=msk)                  ! indexes
      l_node%state%thresholds = vector(tree%train%n_feat, source=null_thr) ! thresholds
    end if

    if (r_node%split) then
      r_node%state%n_obs = p_node%state%best(1)%n_r                        ! size
      r_node%state%idx = pack(p_node%state%idx, mask=.not. msk)            ! indexes
      r_node%state%thresholds = vector(tree%train%n_feat, source=null_thr) ! thresholds
    end if

    call update_state_thresholds( &
      p_node, tree%ctrl%min_obs, tree%train%n_feat, tree%train%x(p_node%state%idx, :), &
      tree%train%na(p_node%state%idx, :), l_node, r_node, update)

    ! If no thresholds were found, deallocate the old node state to free memory
    fail = .not. (update(1) .or. update(2))
    if (fail) then
      call update_node_to_terminal(p_node)
      call update_node_to_terminal(l_node)
      call update_node_to_terminal(r_node)
      return
    end if

    ! At least one child has valid thresholds, update the NA masks.
    ! Copy the important features from the parent node (if needed)
    if (update(1)) then
      l_node%state%n_imp = p_node%state%best(1)%n_imp
      l_node%state%imp_feat = p_node%state%best(1)%imp_feat
    else
      call update_node_to_terminal(l_node)
    end if

    if (update(2)) then
      r_node%state%n_imp = p_node%state%best(1)%n_imp
      r_node%state%imp_feat = p_node%state%best(1)%imp_feat
    else
      call update_node_to_terminal(r_node)
    end if

    ! Early exit: if no missing values, reset the old node state to free memory
    if (.not. tree%train%any_na) then
      call update_node_to_terminal(p_node)
      return
    end if

    ! Check if the p_zero matriz changed for the child nodes
    ! P(father) <= eps implies P(child) <= eps
    n_tn = tree%net%n_tn
    j = findloc(tree%net%tn_id == nodeid, value=.true., dim=1)
    update_zero = tree%net%n_zero(n_tn - 1:n_tn) > tree%net%n_zero(j)
    n = size(p_node%state%na_info, 2)

    ! Copy the masks or perform the update (if needed)
    if (update(1)) then
      l_node%state%na_info = p_node%state%best(1)%na_info
      if (update_zero(1)) then
        do j = 1, n
          where (tree%net%p_zero(:, n_tn - 1))
            l_node%state%na_info(:, j) = .false.
          end where
        end do
      end if
    end if
    if (update(2)) then
      r_node%state%na_info = p_node%state%best(1)%na_info
      if (update_zero(2)) then
        do j = 1, n
          where (tree%net%p_zero(:, n_tn))
            r_node%state%na_info(:, j) = .false.
          end where
        end do
      end if
    end if

    ! Reset the old node state to free memory
    call update_node_to_terminal(p_node)
  end subroutine update_node_state

  pure subroutine update_best_list(n_cand, best, try, fail)
    !-----------------------------------------------------------------------------------
    !  Updates the list of best candidate splits for a node.
    !
    ! ARGUMENTS
    !   n_cand [in]    : Integer, number of candidates.
    !   best   [inout] : info_split(:), current best candidates.
    !   try    [in]    : info_split, new candidate.
    !   fail   [inout] : Logical, flag indicating if update failed.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    integer, intent(in) :: n_cand
    type(info_split), intent(inout) :: best(n_cand)
    type(info_split), intent(in) :: try
    logical, intent(inout) :: fail
    integer :: worst_pos

    fail = .true.
    worst_pos = minloc(best%score, dim=1)
    if (best(worst_pos)%score < try%score) then
      best(worst_pos) = try
      fail = .false.
    end if
  end subroutine update_best_list

  pure subroutine update_best_list_global(n_best, best_id, best_score, nodeid, n_try, try, by_node, first, n_cand_found)
    !-----------------------------------------------------------------------------------
    ! Updates the global list of best candidate splits across all nodes.
    !
    ! ARGUMENTS
    !   n_best       [in]    : Integer, maximum number of best candidates.
    !   best_id      [inout] : integer(n_best,2), current best candidates.
    !   best_score   [inout] : real(n_best), current best scores.
    !   nodeId       [in]    : Integer, node id.
    !   n_try        [in]    : Integer, number of new candidates.
    !   try          [in]    : info_split(:), new candidates.
    !   by_node      [in]    : Logical, flag for node-wise update.
    !   first        [inout] : Logical, flag for first update.
    !   n_cand_found [inout] : Integer, number of candidates found.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------------
    integer, intent(in) :: n_try, nodeid, n_best
    integer, intent(inout) :: best_id(n_best, 2)
    real(dp), intent(inout) :: best_score(n_best)
    type(info_split), intent(in) :: try(n_try)
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
      worst_pos = minloc(best_score, dim=1)
      if (try(i)%score > best_score(worst_pos)) then
        best_score(worst_pos) = try(i)%score
        best_id(worst_pos, 1) = nodeid
        best_id(worst_pos, 2) = i
        n_cand_found = n_cand_found + 1
      end if
    end do
  end subroutine update_best_list_global

  pure subroutine get_idx_prob(idx_p, n_obs, n_info, na_info, n_prob, idx, fill, f, update_f)
    !---------------------------------------------------------------------------------
    ! Updates the index information for probability computation in the current node.
    ! This subroutine determines which observations should be used to compute the
    ! probability matrix `P` for the left and right child nodes, based on the
    ! missingness of features and the chosen fill strategy.
    !
    ! ARGUMENTS
    !   idx_p    [inout] : idx_prob, information on missing values.
    !   n_obs    [in]    : Integer, total number of observations.
    !   na_info  [in]    : Logical, mask for missing values and p > eps
    !   n_prob   [in]    : Integer, number of observations with probabilities > eps.
    !   idx      [in]    : Integer, global indices with p > eps.
    !   fill     [in]    : Integer, fill strategy for missing values.
    !   f        [in]    : Integer, splitting feature
    !   update_f [in]    : loglical, if .true. the splitting feature was just
    !                      added to the imp_feat vector
    !
    ! Added: July/2025
    !---------------------------------------------------------------------------------
    implicit none
    type(idx_prob), intent(inout) :: idx_p
    integer, intent(in) :: n_obs, n_info, n_prob, fill, f
    integer, intent(in) :: idx(n_prob)
    logical, intent(inout) :: na_info(n_obs, n_info)
    logical, intent(in) :: update_f

    ! Set the indices to be used to update the probability matrix P
    select case (fill)
    case (0)
      ! Update the mask for missing values
      if (update_f) na_info(idx, n_info) = na_info(idx, n_info) .or. na_info(idx, f)

      ! Only depends on whether the observation has any missing values or not
      idx_p%n_miss(1) = count(na_info(idx, n_info))
      idx_p%n_miss(2) = 0  ! dummy value
      if (idx_p%n_miss(1) > 0) then
        idx_p%idx_any_na = pack(idx, na_info(idx, n_info))
      end if
      if (idx_p%n_miss(1) < n_prob) then
        idx_p%idx = pack(idx,.not. na_info(idx, n_info))
      end if
    case (1)
      ! Update the mask for missing values
      if (update_f) na_info(idx, n_info) = na_info(idx, n_info) .or. na_info(idx, f)

      ! Depends on both: whether the observation has any missing values or not
      ! and on whether the splitting feature is missing
      idx_p%n_miss(1) = count(na_info(idx, n_info))  ! any missing
      idx_p%n_miss(2) = count(na_info(idx, f))       ! Xf missing
      if (idx_p%n_miss(2) > 0) then
        idx_p%idx_na_f = pack(idx, na_info(idx, f))
      end if
      if (idx_p%n_miss(2) < idx_p%n_miss(1)) then
        idx_p%idx_any_na = pack(idx, na_info(idx, n_info) .and. (.not. na_info(idx, f)))
      end if
      if (idx_p%n_miss(1) < n_prob) then
        idx_p%idx = pack(idx,.not. na_info(idx, n_info))
      end if
    case (2)
      ! Only depends on whether the splitting feature is missing or not
      idx_p%n_miss(1) = 0  ! dummy value
      idx_p%n_miss(2) = count(na_info(idx, f))
      if (idx_p%n_miss(2) > 0) then
        idx_p%idx_na_f = pack(idx, na_info(idx, f))
      end if
      if (idx_p%n_miss(2) < n_prob) then
        idx_p%idx = pack(idx,.not. na_info(idx, f))
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
    !   n_miss  [in]    : Integer, number of observations.
    !   y       [in]    : Real(dp), response values for missing.
    !   y_sum   [inout] : Real(dp), sum for left and right nodes.
    !   idx     [in]    : Integer, index of the missing feature.
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
    real(dp), intent(in) :: y(n_miss)
    real(dp), intent(inout) :: y_sum(2)
    integer, intent(in) :: idx(n_miss), crit
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

  pure subroutine update_prob_no_na(argsd, n_idx, idx, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    implicit none
    !---------------------------------------------------------------------------------------
    ! Computes the probability for the left and right node given a feature and a threshold
    ! when there are no missing values in the important features.
    ! WARNING: This subrotuine implicitly assumes that for idx, all prob > eps.
    !
    ! ARGUMENTS:
    !   argsd     [in] : argsDist, distribution related parameters.
    !   n_idx     [in] : Integer, size of idx.
    !   idx       [in] : Integer, global indices to compute P.
    !   n_obs     [in] : Integer, number of observations.
    !   x_f       [in] : Real(dp), the value for the splitting feature.
    !   threshold [in] : Real(dp), the threshold for the splitting feature.
    !   bounds    [in] : Real(dp), lower and upper bounds for the splitting feature.
    !   sigma_f   [in] : Real(dp), sigma for splitting feature.
    !   prob      [in] : probability for the father node.
    !   p        [out] : probability for the left and right nodes.
    !
    ! Added July/2025
    !---------------------------------------------------------------------------------------
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n_obs, n_idx
    integer, intent(in) :: idx(n_idx)
    real(dp), intent(in) :: x_f(n_obs)
    real(dp), intent(in) :: threshold, bounds(2), sigma_f
    real(dp), intent(in) :: prob(n_obs)
    real(dp), intent(inout) :: p(n_obs, 2)
    real(dp) :: p_temp(n_idx)

    ! Compute P([-Inf, threshold])
    p_temp = pdist(argsd, n_idx, threshold, x_f(idx), sigma_f)

    ! If the lower bound is -Inf, the left node covers (-Inf, threshold]
    ! Otherwise, the left node covers [lower, threshold]
    if (bounds(1) <= neg_inf) then
      p(idx, 1) = p_temp
    else
      p(idx, 1) = p_temp - pdist(argsd, n_idx, bounds(1), x_f(idx), sigma_f)
    end if

    ! If the upper bound is Inf, the right node covers (threshold, Inf)
    ! Otherwise, the right node covers [threshold, upper]
    if (bounds(2) >= pos_inf) then
      p(idx, 2) = 1.0_dp - p_temp
    else
      p(idx, 2) = pdist(argsd, n_idx, bounds(2), x_f(idx), sigma_f) - p_temp
    end if

    ! Standardize the probabilities so that they sum to P(Parent)
    p_temp = prob(idx) / sum(p(idx, :), dim=2)
    p(idx, 1) = p(idx, 1) * p_temp
    p(idx, 2) = p(idx, 2) * p_temp
  end subroutine update_prob_no_na

  pure subroutine update_prob_na_fill0(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    !----------------------------
    ! Helper: fill = 0
    ! Added July/2025
    !----------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    type(idx_prob), intent(in) :: idx_p
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: x_f(n_obs), sigma_f, threshold
    real(dp), intent(in) :: bounds(2), prob(n_obs)
    real(dp), intent(inout) :: p(n_obs, 2)

    if (idx_p%n_miss(1) > 0) then
      p(idx_p%idx_any_na, 1) = prob(idx_p%idx_any_na) / 2.0_dp
      p(idx_p%idx_any_na, 2) = prob(idx_p%idx_any_na) / 2.0_dp
      if (idx_p%n_miss(1) == n_obs) return
    end if

    ! Compute the probability for complete cases
    call update_prob_no_na(argsd, size(idx_p%idx), idx_p%idx, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
  end subroutine update_prob_na_fill0

  pure subroutine update_prob_na_fill1(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    !----------------------------
    ! Helper: fill = 1
    ! Added July/2025
    !----------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    type(idx_prob), intent(in) :: idx_p
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: x_f(n_obs), sigma_f, threshold
    real(dp), intent(in) :: bounds(2), prob(n_obs)
    real(dp), intent(inout) :: p(n_obs, 2)
    integer :: i, n_na

    ! CASE 1: Missing values for X_f
    if (idx_p%n_miss(2) > 0) then
      p(idx_p%idx_na_f, 1) = prob(idx_p%idx_na_f) / 2.0_dp
      p(idx_p%idx_na_f, 2) = prob(idx_p%idx_na_f) / 2.0_dp
      if (idx_p%n_miss(2) == n_obs) return
    end if

    ! CASE 2: X_f observed but missing values in X
    if (idx_p%n_miss(2) < idx_p%n_miss(1)) then
      n_na = idx_p%n_miss(1) - idx_p%n_miss(2)
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
      if (idx_p%n_miss(1) == n_obs) return
    end if

    ! Compute the probability for complete cases
    call update_prob_no_na(argsd, size(idx_p%idx), idx_p%idx, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
  end subroutine update_prob_na_fill1

  pure subroutine update_prob_na_fill2(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    !----------------------------
    ! Helper: fill = 2
    ! Added July/2025
    !----------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    type(idx_prob), intent(in) :: idx_p
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: x_f(n_obs), sigma_f, threshold
    real(dp), intent(in) :: bounds(2), prob(n_obs)
    real(dp), intent(inout) :: p(n_obs, 2)

    ! Compute the indexes for missing values using na_f
    if (idx_p%n_miss(2) > 0) then
      p(idx_p%idx_na_f, 1) = prob(idx_p%idx_na_f) / 2.0_dp
      p(idx_p%idx_na_f, 2) = prob(idx_p%idx_na_f) / 2.0_dp
      if (idx_p%n_miss(2) == n_obs) return
    end if
    ! Compute the probability for complete cases
    call update_prob_no_na(argsd, size(idx_p%idx), idx_p%idx, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
  end subroutine update_prob_na_fill2

  pure subroutine update_prob(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, fill, prob, p)
    !---------------------------------------------------------------------------------------
    ! Computes the probability for the left and right node given a feature and a threshold.
    ! Calls the specific function based on fill type.
    ! WARNING: This subroutine implicitly assumes that for idx, all prob > eps.
    !
    ! ARGUMENTS
    !   argsd    [in] : argsDist, distribution related parameters.
    !   idx_p    [in] : type(idx_prob), information on missing values.
    !   n_obs    [in] : Integer, number of observations.
    !   x_f      [in] : Real(dp), the value for the splitting feature.
    !   na_f     [in] : Logical, missing mask for the current feature.
    !   bounds   [in] : Real(dp), lower and upper bounds for the splitting feature.
    !   sigma_f  [in] : Real(dp), sigma for splitting feature.
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
    integer, intent(in) :: n_obs, fill
    real(dp), intent(in) :: x_f(n_obs), sigma_f, threshold
    real(dp), intent(in) :: bounds(2)
    real(dp), intent(in) :: prob(n_obs)
    real(dp), intent(inout) :: p(n_obs, 2)

    select case (fill)
    case (0)
      ! fill = 0: Assing uniform probability for both nodes if any observation in Xi is missing
      call update_prob_na_fill0(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    case (1)
      ! fill = 1: Assigns uniform probability for both nodes if x_f is missing
      !           Assign 0/1 weights when X_f is not missing
      call update_prob_na_fill1(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    case (2)
      ! fill = 2: Assigns probability based weights.
      call update_prob_na_fill2(argsd, idx_p, n_obs, x_f, threshold, bounds, sigma_f, prob, p)
    end select
  end subroutine update_prob

  pure subroutine search_from_center(state, train, ctrl, argsd, try, best, helper, bounds, idx_p, prob)
    !---------------------------------------------------------------------------------
    ! Perform the search for the best split from the middle to the end of the vector.
    ! The direction of the search depends on the 'direc' argument.
    !
    ! ARGUMENTS
    !   state   [inout] :: node_state, node state to update.
    !   train   [in]    :: tree_data, input data.
    !   ctrl    [in]    :: tree_ctrl, control parameters.
    !   argsd   [in]    :: argsdist, distribution related parameters
    !   try     [inout] :: info_split, information for the split candidate
    !   best    [inout] :: info_split, information on the best splits
    !   helper  [inout] :: h_data, helper variables.
    !   bounds  [in]    :: Real(dp), bounds for current feature
    !   idx_p   [in]    :: idx_prob, information on missing values
    !   prob    [in]    :: Real(dp), vector of probabilites for the father node
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
    !--------------------------------------------------------------------------------
    implicit none
    type(node_state), intent(inout) :: state
    type(tree_data), intent(in) :: train
    type(tree_ctrl), intent(in) :: ctrl
    type(argsdist), intent(in) :: argsd
    type(info_split), intent(inout) :: try
    type(info_split), intent(inout) :: best(ctrl%n_cand)
    type(h_data), intent(inout) :: helper
    type(idx_prob), intent(in) :: idx_p
    real(dp), intent(in) :: bounds(2)
    real(dp), intent(in) :: prob(train%n_obs)
    integer :: i, ii
    logical :: fail
    real(dp) :: y_sum(2), perc_comp(2)

    i = helper%n_start
    ! stop if gets to any end of the vector
    loop: do
      ! increment i based on the search direction (left or right)
      i = i + helper%direction
      if (i < 1 .or. i > state%thresholds(try%feat)%nt) return

      ! save the current threshold
      try%thr = state%thresholds(try%feat)%thr(i)

      ! Find the index ii such that x(ii) > threshold
      ! ii > i, always. The equality ii = i + 1 only holds if x contains only
      ! unique values that are not close to each other.
      ii = i
      do while (helper%x(ii + 1) < try%thr)
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
      call update_prob(argsd, idx_p, train%n_obs, train%x(:, try%feat), try%thr, &
                       bounds, train%sigma(try%feat), train%fill, prob, try%p)
      perc_comp = count(try%p > ctrl%min_prob, dim=1) / real(train%n_obs, dp)

      ! When the first split fails (min_obs or perc_comp), any split from here will also fail
      if (.not. minval(perc_comp) > ctrl%min_prop) return

      ! If a split can be made, save the regions for the new nodes
      try%regid(helper%idx_c(1:ii)) = left_id                  ! left node
      try%regid(helper%idx_c(ii + 1:helper%n_comp)) = right_id ! right node

      ! STEP2: update left and right sum using missing indexes
      ! Assign Xmiss to some reg using a proxy
      if (helper%n_miss_f > 0) then
        call assign_missing(try, helper%n_miss_f, train%y(state%idx(helper%idx_m)), y_sum, helper%idx_m, ctrl%crit)
      end if

      ! compute the proxy score and update the best list
      try%score = y_sum(1)**2 / try%n_l + y_sum(2)**2 / try%n_r

      if (state%n_cand_found < ctrl%n_cand) then
        best(state%n_cand_found + 1) = try
        state%n_cand_found = state%n_cand_found + 1
        cycle
      end if

      call update_best_list(ctrl%n_cand, best, try, fail)
      if (.not. fail) state%n_cand_found = state%n_cand_found + 1
    end do loop
  end subroutine search_from_center

  subroutine find_node_splits(state, train, ctrl, bounds, prob, n_zero, p_zero, argsd)
    !---------------------------------------------------------------------------------------------
    ! Finds the best candidate splits for a node based on a proxy improvement measure.
    !
    ! ARGUMENTS
    !   state  [inout] : node_state, node state to update.
    !   train  [in]    : tree_data, input data.
    !   ctrl   [in]    : tree_ctrl, control parameters.
    !   bounds [in]    : Real(dp), bounds for current region
    !   prob   [in]    : Real(dp), vector of probabilites for the candidate node
    !   n_zero [in]    : Integer, number of prob <= eps
    !   p_zero [in]    : Logical, mask for prob <= eps
    !   argsd  [in]    : argsdist, distribution related parameters
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
    !---------------------------------------------------------------------------------------------
    implicit none
    type(node_state), intent(inout) :: state
    type(tree_data), intent(in) :: train
    type(tree_ctrl), intent(in) :: ctrl
    real(dp), intent(in) :: bounds(train%n_feat, 2)
    real(dp), intent(in) :: prob(train%n_obs)
    integer, intent(in) :: n_zero
    logical, intent(in) :: p_zero(train%n_obs)
    type(argsdist), intent(in) :: argsd
    ! Local Variables
    integer :: j, f, n_obs_node, n_prob
    integer :: idx_node(state%n_obs)
    integer, allocatable :: n_miss(:), idx(:)
    type(info_split), allocatable :: best(:)
    type(info_split) :: try
    logical, allocatable :: not_na_f(:)
    logical :: update_f
    type(h_data) :: helper
    type(idx_prob) :: idx_p, null_idx_p

    ! Number of observations and their indices in the current node
    ! Each child node will have more than min_obs observations (previously computed)
    n_obs_node = state%n_obs
    idx_node = state%idx  ! global indexes
    allocate (n_miss(train%n_feat))  ! to avoid error when compiling for mac
    n_miss = count(train%na(idx_node, :), dim=1) ! number of missing cases in the node

    ! Initialize candidate search for the current node
    best = vector(ctrl%n_cand, try)
    state%n_cand_found = 0
    best%score = 0.0_dp
    try%regid = vector(n_obs_node, 0)
    try%p = matrix(train%n_obs, 2, 0.0_dp)

    ! To preserve probability, if the father's node is near zero,
    ! all the probability is kept in the left node.
    if (n_zero > 0) then
      where (p_zero)
        try%p(:, 1) = prob
        try%p(:, 2) = 0.0_dp
      end where
      ! Global: pack non zero indexes (need update to identify NA entries)
      idx = pack([(j, j=1, train%n_obs)], mask=.not. p_zero)
    else
      ! Global: Non zero indexes (need update to identify NA entries)
      idx = [(j, j=1, train%n_obs)]
    end if
    n_prob = size(idx)

    ! Loop over all features to find best splits
    do f = 1, train%n_feat

      ! Skip feature if it has no valid thresholds (pre-calculated)
      ! if nt > 0, then n_obs > min_obs for each child node (pre-calculated)
      if (state%thresholds(f)%nt == 0) cycle

      ! Pre-process complete cases:
      ! Step 1: filter de data using the global indexes
      !         Here we need all complete cases, not only P > eps.
      not_na_f = .not. train%na(idx_node, f)
      helper%n_comp = count(not_na_f)  ! always > 0
      helper%idx_c = pack(idx_node, mask=not_na_f)
      helper%x = train%x(helper%idx_c, f)
      helper%y_cs = train%y(helper%idx_c)
      helper%n_miss_f = n_miss(f)
      if (helper%n_miss_f > 0) then
        ! position of missing values inside the node block
        ! will be used to assing a region to the missing values
        helper%idx_m = pack([(j, j=1, state%n_obs)], mask=train%na(state%idx, f))
      end if

      ! Step 2: save the local indexes, sort and create the cumulative sum
      helper%idx_c = pack([(j, j=1, n_obs_node)], mask=not_na_f)
      call sort_xy(helper%n_comp, helper%x, helper%y_cs, helper%idx_c)
      do j = 2, helper%n_comp
        helper%y_cs(j) = helper%y_cs(j - 1) + helper%y_cs(j)
      end do

      ! Set current feature for the candidate split
      try%feat = f
      try%n_imp = state%n_imp
      try%imp_feat = state%imp_feat

      ! Select only the important features to update P
      ! Only needs update if f was not an important feature
      if (all(state%imp_feat /= f)) then
        call update_imp_features(try%feat, try%n_imp, try%imp_feat)
        ! NA mask now must also include the possibility of f be missing)
        update_f = (train%n_miss(f) > 0)
      else
        update_f = .false.
      end if

      ! Initialize idx_p to a null state to prevent uninitialized use
      idx_p = null_idx_p

      ! Copy the mask from the parent node (previously updated to reflect NA and P > eps)
      if (train%any_na) then
        try%na_info = state%na_info
        ! set the indexes in idx_p according to the fill type
        call get_idx_prob(idx_p, train%n_obs, size(try%na_info, 2), try%na_info, n_prob, idx, train%fill, f, update_f)
      else
        ! If there are no NA's, no need to save other indexes
        idx_p%n_miss = 0
        idx_p%idx = idx
      end if

      ! Loop in thresholds: starts from the center and goes to the ends (left first)
      ! Stop going in the current direction if a split leads to a region with probability too small

      if (state%thresholds(f)%nt > 1) then  ! More likely
        ! go left (nt/2 to 1)
        helper%n_start = state%thresholds(f)%nt / 2 + 1
        helper%direction = -1
        call search_from_center(state, train, ctrl, argsd, try, best, helper, bounds(f, :), idx_p, prob)
        ! go right (nt/2 + 1 to nt)
        helper%n_start = state%thresholds(f)%nt / 2
        helper%direction = 1
        call search_from_center(state, train, ctrl, argsd, try, best, helper, bounds(f, :), idx_p, prob)
      else
        ! go right (only 1 to test)
        helper%n_start = 0
        helper%direction = 1
        call search_from_center(state, train, ctrl, argsd, try, best, helper, bounds(f, :), idx_p, prob)
      end if
    end do  ! end of loop in features

    ! check if any split can be attempted.
    ! - if n_cand_found == 0, then no split is possible and the node becomes terminal
    if (state%n_cand_found == 0) return

    ! set the best candidates and update the node updating status
    state%n_cand_found = min(state%n_cand_found, ctrl%n_cand)
    state%best = best(1:state%n_cand_found)
    state%update = .false.
  end subroutine find_node_splits

  subroutine update_mse(tree, nodeid, splitid, p, gammahat, yhat, mse)
    !-------------------------------------------------------------------------------------------------
    ! Efficiently updates the coefficients gammahat, the predictions yhat, and the proxy mean
    ! squared error (proxy_mse) for a given split in the tree structure.
    !
    ! ARGUMENTS
    !   tree     [in]    : tree_model, current tree structure.
    !   nodeid   [in]    : Integer, index of the node to split.
    !   splitid  [in]    : Integer, index of the split to apply.
    !   p        [inout] : Real(dp), workspace for the probability matrix (n_obs x (n_tn + 1)).
    !   gammahat [out]   : Real(dp), the updated coefficients.
    !   yhat     [out]   : Real(dp), the updated predictions.
    !   mse      [out]   : Real(dp), the updated proxy mse.
    !
    ! Added: July/2025
    !-------------------------------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(in) :: tree
    integer, intent(in) :: nodeid, splitid
    real(dp), intent(inout) :: p(tree%train%n_obs, tree%net%n_tn + 1)
    real(dp), allocatable, intent(out) :: gammahat(:)
    real(dp), intent(out) :: yhat(tree%train%n_obs)
    real(dp), intent(out) :: mse
    integer :: pos

    ! Finds which column of P corresponds to the node being split
    if (tree%net%n_tn > 1) then
      pos = findloc(tree%net%tn_id == nodeid, value=.true., dim=1)
    else
      pos = 1
    end if

    ! Fill P with old values, except the ones from current region (jth = pos)
    !  - Previous P has n_tn columns.
    !    The columns are indexed by terminal nodes in increasing order
    !  - In the new P, the first n_tn - 1 columns correspond to the old terminal nodes
    !    (which remain the same) and the last two columns correspond to the new regions.
    !    After the jth column, the columns will be shifted one position to the left.
    !  - If P has only one column, there are no columns to shift.
    if (tree%net%n_tn > 1) then
      p(:, 1:pos - 1) = tree%net%p(:, 1:pos - 1)
      p(:, pos:tree%net%n_tn - 1) = tree%net%p(:, pos + 1:tree%net%n_tn)
    end if
    p(:, tree%net%n_tn:tree%net%n_tn + 1) = tree%net%nodes(nodeid)%state%best(splitid)%p

    ! Update gammahat and the proxy mse value
    gammahat = lsquare(tree%train%n_obs, tree%net%n_tn + 1, p, tree%train%y)
    yhat = matmul(p, gammahat)
    mse = sum((tree%train%y - yhat)**2)
  end subroutine update_mse

  pure function expand_nodes(n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a tree_node vector of length n + 2 filled with node at the first
    ! positions and unitialized at the end.
    !
    ! ARGUMENTS
    !   n      [in] : Integer, length of the vector.
    !   source [in] : Type(tree_node), value to fill the vector.
    !
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    type(tree_node), intent(in) :: source(n)
    type(tree_node) :: f(n + 2)
    f(1:n) = source
  end function expand_nodes

  pure function expand_tn(n, j, source1, source2) result(f)
    !---------------------------------------------------------------------------
    ! Creates an integer vector of length n + 1 filled with source1, except the
    ! jth position and source2 at the end.
    !
    ! ARGUMENTS
    !   n       [in] : Integer, length of the vector.
    !   j       [in] : Integer, the position to be replaced
    !   source1 [in] : Integer, values to fill the vector.
    !   source2 [in] : Integer, values to fill the vector.
    !
    ! Added July/2025
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, j
    integer, intent(in) :: source1(n), source2(2)
    integer :: f(n + 1)
    ! Reconstruct the array by taking slices before and after the j-th element,
    ! and then appending the new elements from source2.
    f(1:j - 1) = source1(1:j - 1)
    f(j:n - 1) = source1(j + 1:n)
    f(n:n + 1) = source2
  end function expand_tn

  pure function expand_p(m, n, j, source1, source2) result(f)
    !----------------------------------------------------------------------------
    ! Creates a real(dp) matrix of size m x n filled with source1, except the
    ! jth position and source2 at the end.
    !
    ! ARGUMENTS
    !   m       [in] : Integer, number of rows.
    !   n       [in] : Integer, number of columns.
    !   source1 [in] : Real(dp), value to fill the matrix.
    !   source2 [in] : Real(dp), value to fill the matrix.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, j
    real(dp), intent(in) :: source1(m, n)
    real(dp), intent(in) :: source2(m, 2)
    real(dp) :: f(m, n + 1)
    ! Reconstruct the matrix by taking slices before and after the j-th column,
    ! and then appending the new columns from source2.
    f(:, 1:j - 1) = source1(:, 1:j - 1)
    f(:, j:n - 1) = source1(:, j + 1:n)
    f(:, n:n + 1) = source2
  end function expand_p

  pure function expand_pz(m, n, j, source1, source2) result(f)
    !----------------------------------------------------------------------------
    ! Creates a logical matrix of size m x n filled with source1, except the
    ! jth position and the masks for source2 at the end.
    !
    ! ARGUMENTS
    !   m       [in] : Integer, number of rows.
    !   n       [in] : Integer, number of columns.
    !   source1 [in] : Logical, value to fill the matrix.
    !   source2 [in] : real(dp), value to fill the matrix.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, j
    logical, intent(in) :: source1(m, n)
    real(dp), intent(in) :: source2(m, 2)
    logical :: f(m, n + 1)
    ! Reconstruct the matrix by taking slices before and after the j-th column,
    ! and then appending the new columns from source2.
    f(:, 1:j - 1) = source1(:, 1:j - 1)
    f(:, j:n - 1) = source1(:, j + 1:n)
    f(:, n:n + 1) = source2 <= eps
  end function expand_pz

  pure subroutine update_best_tree_net(tree, split_id)
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
    integer :: newnode(2), n_nodes, n_tn, feat, j
    real(dp) :: thr
    type(info_split) :: best(1)
    type(node_state) :: null_state

    nodeid = split_id(1)
    splitid = split_id(2)
    n_nodes = tree%net%n_nodes
    n_tn = tree%net%n_tn
    feat = tree%net%nodes(nodeid)%state%best(splitid)%feat
    thr = tree%net%nodes(nodeid)%state%best(splitid)%thr

    ! Update P, P_zero, n_zero, gamma, yhat and tn_id before expanding the net
    if (n_tn > 1) then
      j = findloc(tree%net%tn_id == nodeid, value=.true., dim=1)
      tree%net%tn_id = expand_tn(n_tn, j, tree%net%tn_id, [n_nodes + 1, n_nodes + 2])
      tree%net%p = expand_p(tree%train%n_obs, tree%net%n_tn, j, tree%net%p, &
                            tree%net%nodes(nodeid)%state%best(splitid)%p)
      tree%net%p_zero = expand_pz(tree%train%n_obs, tree%net%n_tn, j, tree%net%p_zero, &
                                  tree%net%nodes(nodeid)%state%best(splitid)%p)
      tree%net%n_zero = expand_tn(n_tn, j, tree%net%n_zero, count(tree%net%p_zero(:, n_tn:n_tn + 1), dim=1))
    else
      tree%net%tn_id = [2, 3]
      tree%net%p = tree%net%nodes(nodeid)%state%best(splitid)%p
      tree%net%p_zero = tree%net%nodes(nodeid)%state%best(splitid)%p <= eps
      tree%net%n_zero = count(tree%net%p_zero, dim=1)
    end if

    ! expand the net by adding two new nodes
    tree%net%nodes = expand_nodes(n_nodes, tree%net%nodes)
    tree%net%n_nodes = n_nodes + 2
    tree%net%n_tn = n_tn + 1

    ! Update the old node information
    tree%net%nodes(nodeid)%isterminal = 0
    tree%net%nodes(nodeid)%feature = feat
    tree%net%nodes(nodeid)%threshold = thr
    tree%net%nodes(nodeid)%split = .false.

    ! Update the new nodes basic information
    newnode = [1, 2] + n_nodes
    tree%net%nodes(newnode)%id = [left_id, right_id] + n_nodes
    tree%net%nodes(newnode)%isterminal = 1
    tree%net%nodes(newnode)%fathernode = nodeid
    tree%net%nodes(newnode)%depth = tree%net%nodes(nodeid)%depth + 1
    tree%net%nodes(newnode)%feature = feat
    tree%net%nodes(newnode)%threshold = thr
    tree%net%nodes(newnode)%split = .true.

    ! copy and update the bounds for the new nodes
    !  - update upper bound for the left child
    !  - update lower bound for the right child
    tree%net%nodes(newnode(1))%bounds = tree%net%nodes(nodeid)%bounds
    tree%net%nodes(newnode(2))%bounds = tree%net%nodes(nodeid)%bounds
    tree%net%nodes(newnode(1))%bounds(feat, 2) = thr
    tree%net%nodes(newnode(2))%bounds(feat, 1) = thr

    ! save only the best candidate information for the parent node
    tree%net%nodes(nodeid)%state%n_cand_found = 1
    best(1) = tree%net%nodes(nodeid)%state%best(splitid)
    tree%net%nodes(nodeid)%state%best = best

    ! update the region variable in the net
    tree%net%region(tree%net%nodes(nodeid)%state%idx) = best(1)%regid + n_nodes

    ! Initializing the state variable. New nodes always need update
    ! (other variables will be set after checking basic stopping criterias, if needed)
    tree%net%nodes(newnode(1))%state = null_state
    tree%net%nodes(newnode(2))%state = null_state
    tree%net%nodes(newnode(1))%state%update = .true.
    tree%net%nodes(newnode(2))%state%update = .true.

    ! reset counter
    tree%net%n_cand_found = 0
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
    !   best_id      [in]    : Integer(:,:), IDs of the best candidates.
    !   n_cand_found [in]    : Integer, number of candidates to check.
    !
    ! Added July/2025
    !-----------------------------------------------------------------------------
    implicit none
    type(tree_model), intent(inout) :: tree
    integer, intent(in) :: n_cand_found
    integer, intent(in) :: best_id(n_cand_found, 2)
    real(dp) :: mse_temp, mse_best
    real(dp) :: p(tree%train%n_obs, tree%net%n_tn + 1)
    real(dp), allocatable :: gamma_temp(:), gamma_best(:)
    real(dp) :: yhat_temp(tree%train%n_obs), yhat_best(tree%train%n_obs)
    integer :: k, nodeid, splitid, split_id(2)

    mse_best = tree%net%mse ! sum of squares of residuals
    split_id = -1

    do k = 1, n_cand_found
      nodeid = best_id(k, 1)
      splitid = best_id(k, 2)

      ! Compute the mse for current candidate and check if the new split improves the mse
      ! Recompute gammahat and yhat for the best tree at the end
      call update_mse(tree, nodeid, splitid, p, gamma_temp, yhat_temp, mse_temp)

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
      tree%net%nodes%split = .false.
      return
    end if

    ! If a split was found, check if the split is worth it using cp criterion.
    ! If the split is not worth it, disables further splitting and returns.
    if (mse_best > (1 - tree%ctrl%cp) * tree%net%mse) then
      tree%net%nodes%split = .false.
      return
    end if

    ! full update the tree using the best candidate
    tree%net%gammahat = gamma_best
    tree%net%yhat = yhat_best
    tree%net%mse = mse_best
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

    ! Main loop to create divisions.
    !  - Uses a stopping criteria based on the number of terminal nodes (max_tn)
    grow_loop: do while ((tree%net%n_tn < tree%ctrl%max_tn))

      ! No need to check when n_tn = 1 (always ok)
      ! If n_tn > 1 update do_split and n_split using basic stopping criteria
      !  - criteria used in this step: depth, min_prop, min_prob and max_d
      !    (only needs to check the last two childs)
      if (tree%net%n_tn == 1) then
        n_split = 1
        do_split = [1]
      else
        ! Updates the list of splittable nodes
        call update_split_candidates(tree%net, tree%ctrl, tree%train%n_obs, n_split, do_split)
        ! If no further splits are possible we are done growing the three
        if (n_split == 0) then
          if (printinfo >= log_detailed) then
            ! Debug: Print node splitting info
            call labelpr("    No more candidates to split", -1)
            call labelpr(" ", -1)
          end if
          return
        end if

        ! If n_tn = 1 thresholds and NA masks were already computed during initialization
        ! Otherwhise, update the thresholds information and NA masks for the nodes created in the last iteration
        ! (to avoid unnecessary updates child nodes were not updated when the parent node was splitted)
        call update_node_state(tree, fail)

        ! if the last two nodes are the only splittable nodes and there are no more
        ! thresholds then the search for a new tree is over (no need to update to terminal)
        if (fail) then
          if (n_split <= 2 .and. do_split(1) >= tree%net%n_nodes - 1) return
        end if
      end if

      ! Atempting a new split.
      !  - If n_tn = 1, finds the best candidates
      !  - if n_tn > 1 loop over the new nodes (if splittable) and update the list of candidates
      if (printinfo >= log_detailed) then
        ! Debug: Print node splitting info
        call labelpr("**************************************************", -1)
        call labelpr("    Attempting a new split ", -1)
        call labelpr("**************************************************", -1)
        call labelpr("Current status", -1)
        call intpr1("    Nodes (n_nodes):", -1, tree%net%n_nodes)
        call intpr1("    Terminal nodes (n_tn):", -1, tree%net%n_tn)
        call intpr("    Candidates to split:", -1, do_split, n_split)
        call labelpr(" ", -1)
      end if

      ! Setting variables used int the search for the best candidates
      ! - by_node: controls how the search is done,
      !     if .true., n_cand candidates for each node
      !     if .false., n_cand candidates globally
      ! - first: controls if best_id is to be initialized or updated
      ! - best_id and best_score: informations about the best candidates

      n_best = tree%ctrl%n_cand          ! at least 1 candidate is needed
      if (tree%ctrl%by_node) n_best = n_best * n_split
      best_id = matrix(n_best, 2, 0)
      best_score = vector(n_best, 0.0_dp)
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
        if (.not. tree%net%nodes(nodeid)%split) cycle

        ! Check if the the node needs update
        ! If update = .false. then the best candidates for this node remain the same
        if (.not. tree%net%nodes(nodeid)%state%update) then
          ! Update the best candidates list
          call update_best_list_global( &
            n_best, best_id, best_score, nodeid, tree%net%nodes(nodeid)%state%n_cand_found, &
            tree%net%nodes(nodeid)%state%best, tree%ctrl%by_node, first, n_cand_found)
          cycle
        end if

        if (printinfo >= log_detailed) then
          ! Debug: Print node splitting info
          call intpr1("    Searching candidates for node Id:", -1, nodeid)
          call intpr1("    Node observations (n_obs):", -1, tree%net%nodes(nodeid)%state%n_obs)
          if (printinfo >= log_debug_deep) then
            call intpr("    Indexes", -1, tree%net%nodes(nodeid)%state%idx, tree%net%nodes(nodeid)%state%n_obs)
          end if
          call labelpr(" ", -1)
        end if

        ! loop over all variables to find the best n_cand for the current node
        col_id = findloc(tree%net%tn_id == nodeid, value=.true., dim=1)
        call find_node_splits(tree%net%nodes(nodeid)%state, tree%train, tree%ctrl, &
                              tree%net%nodes(nodeid)%bounds, tree%net%p(:, col_id), &
                              tree%net%n_zero(col_id), tree%net%p_zero(:, col_id), tree%dist)

        ! disabling further split attempts for nodes that do not have enough data.
        if (tree%net%nodes(nodeid)%state%n_cand_found == 0) then
          call update_node_to_terminal(tree%net%nodes(nodeid))
          cycle
        end if

        ! if at least one candidate was found, update the best candidates list
        call update_best_list_global( &
          n_best, best_id, best_score, nodeid, tree%net%nodes(nodeid)%state%n_cand_found, &
          tree%net%nodes(nodeid)%state%best, tree%ctrl%by_node, first, n_cand_found)
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
  pure function prob_rgivenx_no_na(argsd, n_obs, x, bounds, sigma) result(prob)
    !-------------------------------------------------------------------------------------
    ! Compute the probability P(R | X) = Psi(X, R, sigma) when there are no missing values
    ! and only one feature. Uses the vectorized form of pdist
    !
    ! ARGUMENTS
    !   argsd   [in] : argsDist, distribution related parameters.
    !   n_obs   [in] : Integer, number of observations.
    !   X       [in] : Real(dp), vector of feature values (n_obs).
    !   bounds  [in] : Real(dp), region bounds.
    !   sigma   [in] : Real(dp), standard deviations for each feature.
    !
    ! Added July/2025
    !-------------------------------------------------------------------------------------
    implicit none
    type(argsdist), intent(in) :: argsd
    integer, intent(in) :: n_obs
    real(dp), intent(in) :: x(n_obs)
    real(dp), intent(in) :: bounds(2)
    real(dp), intent(in) :: sigma
    real(dp) :: prob(n_obs)

    ! The conditional checks on bounds are structured to handle the most common case (finite bounds) first.
    if (bounds(1) > neg_inf .and. bounds(2) < pos_inf) then
      ! Case 1: R = [lower, upper]. (Most likely)
      prob = pdist(argsd, n_obs, bounds(2), x, sigma) - pdist(argsd, n_obs, bounds(1), x, sigma)
    else if (bounds(1) <= neg_inf) then
      ! Case 2: Lower bound is -inf
      if (bounds(2) < pos_inf) then
        ! Subcase 2.1: R = (-inf, upper].
        prob = pdist(argsd, n_obs, bounds(2), x, sigma)
      else
        ! Subcase 2.2: R = (-inf, +inf). (Very unlikely)
        prob = 1.0_dp
      end if
    else
      ! Case 3: R = [lower, +inf).
      prob = 1.0_dp - pdist(argsd, n_obs, bounds(1), x, sigma)
    end if
  end function prob_rgivenx_no_na

  pure subroutine prob_rgivenx_sf(n_obs, x, na, n_tn, depth, bounds, p, sigma, argsd)
    implicit none
    integer, intent(in) :: n_obs, n_tn
    integer, intent(in) :: depth(n_tn)
    real(dp), intent(in) :: x(n_obs), bounds(n_tn, 2), sigma
    logical, intent(in) :: na(n_obs)
    real(dp), intent(inout) :: p(n_obs, n_tn)
    type(argsdist), intent(in) :: argsd
    integer :: i, k, n_comp
    logical :: is_complete
    integer, allocatable :: idx(:)

    ! check if the data has missing values
    is_complete = .not. any(na)

    ! DIRECT APPROACH: Optimization for complete data and single feature
    if (is_complete) then
      do k = 1, n_tn
        p(:, k) = prob_rgivenx_no_na(argsd, n_obs, x, bounds(k, :), sigma)
      end do
      return
    end if

    ! (a) Process all entries that do not have missing values.
    idx = pack([(i, i=1, n_obs)], mask=.not. na)
    n_comp = size(idx)
    if (n_comp > 0) then
      ! DIRECT APPROACH: Optimization for complete data and single feature.
      ! The maximum number of calls to pdist is 2 * n_tn * n_obs
      do k = 1, n_tn
        p(idx, k) = prob_rgivenx_no_na(argsd, n_comp, x(idx), bounds(k, :), sigma)
      end do
    end if

    ! (b) If no entries with missing values remain, exit.
    if (n_comp == 0) then
      idx = [(i, i=1, n_obs)]
    else
      idx = pack([(i, i=1, n_obs)], mask=na)
    end if

    ! (c) Process missing entries hierarchically. Since there exist only one feature,
    !     P(Left) = P(Right) = P(Father)/2 = 1/2^depth
    !     This is calculated directly for each terminal node.
    do k = 1, n_tn
      if (depth(k) < max_prob_depth) then
        p(idx, k) = 1.0_dp / 2**depth(k)
      else
        p(idx, k) = 0.0_dp
      end if
    end do
  end subroutine prob_rgivenx_sf

  pure function prob_rgivenx(argsd, n_obs, n_feat, x, sigma, fill, n_tn, tn, nodes_info, bounds) result(p)
    !-----------------------------------------------------------------------------
    ! Computes probability matrix P for new data with missing values.
    ! The calculation is done hierarchically to mimic the procedure followed
    ! during the tree building process.
    !
    ! ARGUMENTS
    !   argsd      [in] : argsDist, distribution related parameters.
    !   n_obs      [in] : Integer, number of observations.
    !   n_feat     [in] : Integer, number of features.
    !   X          [in] : Real(dp), data matrix (n_obs x n_feat).
    !   sigma      [in] : Real(dp), vector of sigmas (n_feat).
    !   fill       [in] : Integer, fill strategy for missing values.
    !   n_tn       [in] : Integer, number of terminal nodes.
    !   tn         [in] : Integer, terminal nodes index.
    !   nodes_info [in] : Integer, information on nodes. (2*n_tn - 1, 4)
    !                     (isTerminal, father, depth, feature)
    !   bounds     [in] : Real(dp), bounds for each feature in the nodes.
    !                     (n_feat * (2 * n_tn - 1), 2)
    !   p          [out]: Real(dp), the resulting probability matrix
    !                     of size (n_obs x n_tn).
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
    integer, intent(in) :: n_obs, n_feat, n_tn, fill
    integer, intent(in) :: tn(n_tn), nodes_info(2 * n_tn - 1, 4)
    real(dp), intent(in) :: x(n_obs, n_feat), sigma(n_feat)
    real(dp), intent(in) :: bounds(n_feat * (2 * n_tn - 1), 2)
    real(dp) :: p(n_obs, n_tn)

    ! Local variables
    integer :: i, k, n_nodes, father_id, f, n_miss(n_feat)
    integer :: n_imp, n_prob, n_zero, n_info
    integer, allocatable :: imp_feat(:), idx(:)
    logical :: na(n_obs, n_feat), is_complete, p_zero(n_obs)
    logical, allocatable :: na_info(:, :)
    real(dp) :: p_matrix(n_obs, 2 * n_tn - 1)
    real(dp) :: bd(n_feat, 2), threshold, bd1(n_tn, 2)
    type(idx_prob) :: idx_p, null_idx_p

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
      idx = [(i, i=1, n_feat)]
      ! Features used to build the tree
      imp_feat = pack(idx, mask=[(any(idx(i) == nodes_info(:, 4)), i=1, n_feat)])
      n_imp = size(imp_feat)
      ! mask for complete dataset
      is_complete = all(n_miss(imp_feat) == 0)
    else
      n_imp = 1
      imp_feat = [1]
    end if

    ! Case 1: Single important feature.
    !  - find the indexes for complete data and use the direct approach.
    !  - find the indexes for missing data and use the hierarchical approach (if any)
    if (n_imp == 1) then
      f = imp_feat(1)
      do i = 1, n_tn
        bd1(i, :) = bounds((tn(i) - 1) * n_feat + f, :)
      end do
      call prob_rgivenx_sf(n_obs, x(:, f), na(:, f), n_tn, nodes_info(tn, 3), bd1, p, sigma(f), argsd)
      return
    end if

    ! Case 2: Multiple features. Process all entries hierarchically.
    ! HIERARCHICAL APPROACH:  The maximum number of calls to pdist is 3 * (n_tn - 1) * n_obs.
    ! Loop through sibling pairs to compute probabilities.

    ! Allocate na_info once before the loop to avoid repeated allocations
    if (.not. is_complete) then
      n_info = merge(n_feat, n_feat + 1, fill == 2)
      na_info = matrix(n_obs, n_info, .false.)
    end if

    p_matrix(:, 1) = 1.0_dp  ! root node
    do k = 2, n_nodes - 1, 2
      ! find the current father node id and splitting feature
      father_id = nodes_info(k, 2)
      f = nodes_info(father_id, 4)

      ! To preserve probability, if the father's node is near zero,
      ! all the probability is kept in the left node.
      p_zero = (p_matrix(:, father_id) <= eps)
      n_zero = count(p_zero)
      if (n_zero > 0) then
        where (p_zero)
          p_matrix(:, k) = p_matrix(:, father_id)
        end where
        ! pack non zero indexes (needs update to remove NA entries)
        idx = pack([(i, i=1, n_obs)], mask=.not. p_zero)
      else
        ! pack non zero indexes (needs update to remove NA entries)
        idx = [(i, i=1, n_obs)]
      end if
      n_prob = size(idx)

      ! bounds for the father node and splitting threshold
      bd = bounds(n_feat * (father_id - 1) + 1:n_feat * father_id, :)
      threshold = bounds(n_feat * (k - 1) + f, 2)

      ! select important features for the current node using bounds for the
      ! father node and the current feature id.
      call select_imp_features(f, n_feat, bd, imp_feat, n_imp)

      ! Initialize idx_p to a null state to prevent uninitialized use
      idx_p = null_idx_p

      ! if there are no missing values, then only checks the p_zero mask
      if (is_complete) then
        idx_p%n_miss = 0
        idx_p%idx = idx ! p > 0
      else
        ! update the masks to take into account p_zero.
        ! This changes for each node and cannot be set outside the loop
        na_info(:, 1:n_feat) = na
        do i = 1, n_feat
          where (p_zero)
            na_info(:, i) = .false. ! Xf could be missing but P(father) <= eps
          end where
        end do
        ! Set the indexes bases on the masks and the fill type
        select case (fill)
        case (0, 1)
          ! Find the indexes where any X is missing and P > eps
          na_info(:, n_info) = any(na_info(:, imp_feat), dim=2)
          call get_idx_prob(idx_p, n_obs, n_info, na_info, n_prob, idx, fill, f, .false.)
        case (2)
          call get_idx_prob(idx_p, n_obs, n_info, na_info, n_prob, idx, fill, f, .false.)
        end select
      end if

      ! update the probability matrix
      call update_prob(argsd, idx_p, n_obs, x(:, f), threshold, bd(f, :), sigma(f), fill, &
                       p_matrix(:, father_id), p_matrix(:, k:k + 1))
    end do
    p = p_matrix(:, tn)
  end function prob_rgivenx

  pure subroutine predict_tree(tree)
    implicit none
    type(tree_model), intent(inout) :: tree
    integer :: n_feat, i
    integer :: nodes_info(tree%net%n_nodes, 4)
    real(dp) :: bounds(tree%train%n_feat * tree%net%n_nodes, 2)
    real(dp) :: yhat(tree%test%n_obs)

    ! Set the information to pass to the probability function
    n_feat = tree%train%n_feat
    nodes_info(:, 1) = tree%net%nodes%isterminal
    nodes_info(:, 2) = tree%net%nodes%fathernode
    nodes_info(:, 3) = tree%net%nodes%depth
    nodes_info(:, 4) = tree%net%nodes%feature

    ! Find the bounds for each node and important features
    do i = 1, tree%net%n_nodes
      bounds((i - 1) * n_feat + 1:i * n_feat, :) = tree%net%nodes(i)%bounds
    end do

    !---------------------------------------------
    ! filling the matrix P with the values
    !    Psi(x, Rj, sigma), 1 <= j <= n_tn
    ! where Rj is the j-th terminal node (reg)
    !---------------------------------------------
    tree%test%p = prob_rgivenx(tree%dist, tree%test%n_obs, n_feat, tree%test%x, tree%train%sigma, &
                               tree%train%fill, tree%net%n_tn, tree%net%tn_id, nodes_info, bounds)
    yhat = matmul(tree%test%p, tree%net%gammahat)
    tree%test%yhat = yhat
    tree%test%mse = sum((tree%test%yhat - tree%test%y)**2)
  end subroutine predict_tree

  !===============================================================================
  ! STEP 5: Extract the information from the treee object
  !===============================================================================
  pure subroutine return_tree_train(net, n_obs, n_feat, n_train, idx_train, max_tn, n_tn, &
                               nodes_info, thresholds, bounds, p, gammahat, yhat, mse, xregion)
    !----------------------------------------------------------------------------
    ! Returns the final tree structure to R.
    ! This subroutine copies the values from the final tree structure
    ! to the output vectors and matrices.
    !
    ! ARGUMENTS
    !   net        [in]    : tree_net, final tree structure.
    !   n_obs      [in]    : Integer, number of observations.
    !   n_feat     [in]    : Integer, number of features.
    !   n_train    [in]    : Integer, number of training observations.
    !   idx_train  [in]    : Integer(:), indices of training observations.
    !   max_tn     [in]    : Integer, maximum number of terminal nodes.
    !   n_tn       [out]   : Integer, number of terminal nodes in the final tree.
    !   nodes_info [out]   : Integer, matrix with information about the nodes
    !                        in the final tree. (2 * max_tn - 1, 5)
    !   thresholds [inout] : Real(dp), vector of thresholds defining the
    !                         regions in the final tree. (2 * max_tn - 1)
    !   bounds     [inout] : Real(dp), bounds for each feature in the nodes.
    !                        (n_feat * (2 * max_tn - 1), 2)
    !   P          [inout] : Real(dp), probability matrix for the final tree
    !                        (n_obs x max_tn).
    !   gammahat   [inout] : Real(dp), vector of coefficients for the final tree
    !                        (max_tn).
    !   yhat       [inout] : Real(dp), predicted values for Y (n_obs).
    !   mse        [inout] : Real(dp), mean square error for the final tree.
    !   Xregion    [inout] : Integer, matrix with the indices of the features
    !                        used in the nodes (n_obs).
    !----------------------------------------------------------------------------
    implicit none
    type(tree_net), intent(in) :: net
    integer, intent(in) :: n_obs, n_feat, max_tn, n_train
    integer, intent(in) :: idx_train(n_train)
    integer, intent(out) :: nodes_info(2 * max_tn - 1, 5)
    integer, intent(out) :: n_tn
    integer, intent(out) :: xregion(n_obs)
    real(dp), intent(out) :: thresholds(2 * max_tn - 1)
    real(dp), intent(out) :: p(n_obs, max_tn), gammahat(max_tn), yhat(n_obs)
    real(dp), intent(out) :: bounds(n_feat * (2 * max_tn - 1), 2)
    real(dp), intent(out) :: mse
    integer :: n_nodes, k

    ! Updating the number of terminal nodes
    ! number of nodes in the final tree
    n_nodes = net%n_nodes

    ! Updating the number of terminal nodes and P
    n_tn = net%n_tn
    p = 0.0_dp
    do k = 1, n_tn
      p(idx_train, k) = net%p(:, k)
    end do

    ! Updating gammahat
    gammahat = 0.0_dp
    gammahat(1:n_tn) = net%gammahat

    ! predicted values and mse
    yhat(idx_train) = net%yhat
    mse = net%mse / n_train

    ! initializing and updating the nodes_info
    nodes_info = 0
    nodes_info(1:n_nodes, 1) = net%nodes%id
    nodes_info(1:n_nodes, 2) = net%nodes%isterminal
    nodes_info(1:n_nodes, 3) = net%nodes%fathernode
    nodes_info(1:n_nodes, 4) = net%nodes%depth
    nodes_info(1:n_nodes, 5) = net%nodes%feature

    ! initializing and updating the thresholds
    thresholds = 0.0_dp
    thresholds(1:n_nodes) = net%nodes%threshold

    ! initializing and updating the bounds
    bounds = 0.0_dp
    do k = 1, n_nodes
      bounds((k - 1) * n_feat + 1:k * n_feat, 1) = net%nodes(k)%bounds(:, 1)
      bounds((k - 1) * n_feat + 1:k * n_feat, 2) = net%nodes(k)%bounds(:, 2)
    end do

    ! updating Xregion
    xregion(idx_train) = net%region
  end subroutine return_tree_train

  pure subroutine return_tree_test(test, n_obs, n_test, idx_test, n_tn, p, yhat, mse, xregion)
    !----------------------------------------------------------------------------
    ! Returns the final tree structure to R.
    ! This subroutine copies the values from the final tree structure
    ! to the output vectors and matrices.
    !
    ! ARGUMENTS
    !   test     [in]    : tree_tdata, testing data.
    !   n_obs    [in]    : Integer, number of observations.
    !   n_test   [in]    : Integer, number of testing observations.
    !   idx_test [in]    : Integer(:), indices of testing observations.
    !   n_tn     [in]    : Integer, number of terminal nodes in the final tree.
    !   P        [inout] : Real(dp), probability matrix for the final tree
    !                      (n_obs x n_tn).
    !   yhat     [inout] : Real(dp), predicted values for Y (n_obs).
    !   mse      [inout] : Real(dp), mean square error for the final tree.
    !   Xregion  [inout] : Integer, matrix with the indices of the features
    !                        used in the nodes (n_obs).
    !----------------------------------------------------------------------------
    implicit none
    type(tree_tdata), intent(in) :: test
    integer, intent(in) :: n_obs, n_test, n_tn
    integer, intent(in) :: idx_test(n_test)
    integer, intent(inout) :: xregion(n_obs)
    real(dp), intent(inout) :: p(n_obs, n_tn), yhat(n_obs)
    real(dp), intent(out) :: mse
    integer :: i

    ! Updating P
    do i = 1, n_tn
      p(idx_test, i) = test%p(:, i)
    end do

    ! predicted values and mse
    yhat(idx_test) = test%yhat
    mse = test%mse / n_test

    ! For thest sample, there is no region associated to X.
    xregion(idx_test) = 0
  end subroutine return_tree_test
end module prtree_main
