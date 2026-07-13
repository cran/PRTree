! Contains
!  - Module parameters and Derived types
!  - Interfaces to R functions
!  - Initializations subroutines (to allocate vector and matrices)

module prtree_types
  implicit none
  !=======================================
  ! Module parameters
  !=======================================
  ! Double precision kind
  integer, parameter :: dp = kind(1.0d0)

  ! Infinity for double precision
  real(dp), parameter :: pos_inf = huge(1.0_dp)   ! Positive
  real(dp), parameter :: neg_inf = -huge(1.0_dp)  ! Negative
  real(dp), parameter :: eps = epsilon(1.0_dp)    ! Machine epsilon

  ! Threshold for feature importance
  real(dp), parameter :: feature_threshold = 1.0e-7_dp

  ! Maximum depth for direct probability calculation 1/2^depth.
  ! This is derived from machine epsilon to avoid magic numbers. A depth
  ! greater than this results in a probability smaller than epsilon(1.0_dp),
  ! which is numerically negligible. For double precision, this is 53.
  integer, parameter :: max_prob_depth = exponent(1.0_dp / eps) + 1

  ! Logging levels for printInfo
  integer, parameter :: log_quiet = -1              ! No output
  integer, parameter :: log_base = 0                ! Basic output (errors, final results)
  integer, parameter :: log_verbose = 1             ! More detailed info on splits
  integer, parameter :: log_detailed = 2            ! Verbose info on candidates and MSE
  integer, parameter :: log_debug = 3               ! Full debug output (indexes, thresholds)
  integer, parameter :: log_debug_deep = 10         ! Deep debug

  ! ids for left and righ nodes. Do not change this values        !
  integer, parameter :: left_id = 1          ! id for the left node
  integer, parameter :: right_id = 2         ! id for the right node

  ! Column indices for nodes_info matrix
  integer, parameter :: idx_ni_id = 1
  integer, parameter :: idx_ni_terminal = 2
  integer, parameter :: idx_ni_father = 3
  integer, parameter :: idx_ni_depth = 4
  integer, parameter :: idx_ni_feature = 5

  ! The current distributions available
  character(len=8), parameter :: dist1 = "norm"
  character(len=8), parameter :: dist2 = "lnorm"
  character(len=8), parameter :: dist3 = "gamma"
  character(len=8), parameter :: dist4 = "t"
  ! array with the codes for each distribution
  character(len=8), parameter :: current_dist(4) = [dist1, dist2, dist3, dist4]

  !=================================
  ! Derived types
  !=================================
  ! tree_model: the main tree object
  !  = > argsDist: distribution related information
  !  = > tree_data: information on the input data
  !   == > tree_node: individual tree nodes
  !      == > info_split: information on splitting
  !=================================
  type idx_prob
    !---------------------------------------------
    ! helper variable to simplify functions calls
    ! Added July/2025
    ! Modified March/2026
    !---------------------------------------------
    ! indexes with P > eps
    integer :: n_prob = 0
    integer, allocatable :: idx(:)

    ! indexes with Xf miss and P > eps
    integer :: n_miss_f = 0
    integer, allocatable :: idx_na_f(:)

    ! indexes with X miss and P > eps
    integer :: n_miss_any = 0
    integer, allocatable :: idx_any_na(:)
  end type idx_prob

  type info_split
    !----------------------------------------------------------------
    ! information on n_cand best candidates
    ! For each node, use one for each candidate
    ! Added: July/2025
    !----------------------------------------------------------------
    integer :: feat = -1             ! feature
    real(dp) :: thr = 0.0_dp         ! threshold
    integer :: n_l = 0               ! left node size
    integer :: n_r = 0               ! right node size
    real(dp) :: score = 0.0_dp       ! proxy score for the split
    integer :: n_imp = 0             ! number of important features
    integer, allocatable :: imp_feat(:) ! important features
    integer, allocatable :: regid(:) ! reg number after split
  end type info_split

  type info_thr
    !--------------------------------------------------------------
    ! Information on thresholds for a feature
    ! For each node, use one for each feature
    ! Added: July, 2025
    !--------------------------------------------------------------
    integer :: nt = 0                       ! number of thresholds found
    real(dp), allocatable :: thr(:)         ! thresholds found
  end type info_thr

  type tree_node
    !----------------------------------------------------------------------
    ! Node information.
    ! Use one for each node in the tree
    ! Last update: May/2026
    !  - removed the node_state derived type and added the parameters
    !    directly in the tree_node type
    !----------------------------------------------------------------------
    ! Node statistics and information to help splitting the node
    logical :: split = .true.      ! splitting condition
    logical :: update = .true.     ! if the node needs update
    integer :: n_cand_found = 0    ! number of candidates found in this node
    integer :: n_obs_node = 0      ! node size
    integer :: n_imp = 0           ! number of important features for this node

    ! allocatable arrays: update during the split process
    ! - deallocate after splitting to save memory
    integer, allocatable :: idx(:)               ! indexes of X in this node
    integer, allocatable :: imp_feat(:)          ! important features for this node
    type(info_thr), allocatable :: thresholds(:) ! thresholds for each variable
    type(info_split), allocatable :: best(:)     ! splitting candidates
  end type tree_node

  type argsdist
    !-------------------------------------------------------------
    ! Distribution related parameters
    ! Added July/2025
    !-------------------------------------------------------------
    integer :: dist_id = 1        ! code identifying the distribution
    real(dp) :: dist_par = 0.0_dp ! extra parameter (if any)
  end type argsdist

  type tree_model
    !-------------------------------------------------------------
    ! tree object
    ! Added July/2025
    ! Last update: May/2026
    !   - removed the tree_net/ctrl/tree_data/tree_tdata
    !     derived type and added the parameters directly here
    !-------------------------------------------------------------
    integer :: printinfo = -1     ! logging level
    integer :: n_obs = 0          ! number of observations (training + testing)
    integer :: n_feat = 0         ! number of features (covariates)

    !-----------------------------------------------------------------------
    ! Criteria used to build the tree
    !-----------------------------------------------------------------------
    integer :: fill = 2            ! method used to fill P
    integer :: crit = 3            ! criterion to assign missing
    integer :: max_tn = 1          ! Max terminal nodes
    integer :: max_d = 1           ! Max depth
    integer :: min_obs = 1         ! Min observations
    integer :: n_cand = 3          ! Candidates count
    logical :: by_node = .false.   ! find best candidate by node or global
    real(dp) :: min_prop = 0.1_dp  ! Min proportion
    real(dp) :: min_prob = 0.05_dp ! Min probability
    real(dp) :: cp = 0.01_dp       ! Complexity param (reduction in mse)

    !------------------------------------------------------------------
    ! information on the training and testing data
    !------------------------------------------------------------------
    integer :: n_train = 0                                   ! number of observations in the training data
    real(dp), pointer, contiguous :: y_train(:) => null()    ! Response (n_train)
    real(dp), pointer, contiguous :: x_train(:, :) => null() ! covariates (n_train x n_feat)
    real(dp), pointer, contiguous :: yhat_train(:) => null() ! predictions (n_train)
    real(dp), allocatable :: yhat_temp(:)                    ! helper array (n_train)
    logical, allocatable :: na(:, :)                         ! missing flag (n_train x n_feat)
    integer, allocatable :: n_miss(:)                        ! number of missing cases for each feature
    logical :: any_na = .false.                              ! if the data has any NA
    real(dp) :: mse_train = pos_inf                          ! mean squared errors

    !------------------------------------------------------------------
    ! information on the testing data
    !------------------------------------------------------------------
    integer :: n_test = 0                                    ! number of observations in the testing data
    real(dp), pointer, contiguous :: y_test(:) => null()     ! Response (n_test)
    real(dp), pointer, contiguous :: x_test(:, :) => null()  ! covariates (n_test x n_feat)
    real(dp),  pointer, contiguous :: yhat_test(:) => null() ! predictions (n_test)
    real(dp), allocatable :: yhat_ttemp(:)                   ! helper array (n_test)
    real(dp), pointer, contiguous:: p_test(:, :) => null()   ! matrix P for the test data (n_test x n_tn)
    real(dp), allocatable :: p_test_temp(:, :)               ! helper array (n_test x n_tn)
    real(dp) :: mse_test = pos_inf                           ! mean squared errors

    !------------------------------------------------------------------
    ! Net arguments
    ! Tree nodes and the corresponding parameters and predicitions
    !------------------------------------------------------------------
    integer :: n_nodes = 1                   ! number of nodes
    integer :: n_tn = 1                      ! number of terminal nodes
    integer :: n_cand_found = 0              ! number of candidates found
    integer, allocatable :: tn_id(:)         ! id for terminal node (n_tn)

    integer, pointer :: region(:) => null()  ! region for each X (n_train)
    integer, allocatable :: region_temp(:)   ! helper array for regions

    real(dp), pointer :: p(:, :) => null()   ! matrix P (n_train x n_tn)
    real(dp), allocatable :: p_temp(:, :)    ! helper array for P

    real(dp), pointer :: gammahat(:)         ! coefficients (n_tn)
    real(dp), allocatable :: gamma_temp(:)   ! helper array for coefficients

    type(tree_node), allocatable :: nodes(:)       ! nodes object (n_nodes)
    integer, pointer :: nodes_info(:, :) => null() ! nodes information (n_nodes x 5)
    integer, allocatable :: nodes_info_temp(:, :)  ! helper array for nodes information

    real(dp), pointer :: thresholds(:) => null() ! thresholds for each feature (2 * n_tn -1)
    real(dp), allocatable :: thresholds_temp(:)  ! helper array for thresholds

    !------------------------------------------------------------------
    ! Workspace for least squares (DGELSD)
    !------------------------------------------------------------------
    integer :: lsq_n = 0
    integer :: lsq_lwork = 0
    integer :: lsq_liwork = 0

    !------------------------------------------------------------------
    ! Distribution related parameters
    !------------------------------------------------------------------
    real(dp), allocatable :: sigma(:) ! distribution parameter (n_feat)
    type(argsdist) :: dist
  end type tree_model

  type h_data
    !---------------------------------------------
    ! helper variable to simplify functions calls
    ! Added July/2025
    !---------------------------------------------
    integer :: n_comp = 0            ! local: number of complete cases
    integer :: n_miss_f = 0          ! local: number of missing cases
    integer, allocatable :: idx_c(:) ! local: indexes of complete cases
    integer, allocatable :: idx_m(:) ! local: indexes of missing cases
    real(dp), allocatable :: x(:)    ! local: sorted complete xf
    real(dp), allocatable :: y_cs(:) ! local: sorted y
    real(dp), allocatable :: y_m(:)  ! local: y values for missing cases
    integer :: n_start = 0           ! threshold search
    integer :: direction = 1         ! threshold search
  end type h_data

  !=================================
  ! Interfaces to R functions
  !=================================
  interface
    !--------------------------------------------------------------------
    ! Interfaces for external C functions from R's math library.
    !--------------------------------------------------------------------
    ! Added: July/2025
    ! Modified: March/2026
    !  - transoformed to vectorized versions to avoid loops in Fortran
    !    and speed up the code
    subroutine pnorm_v(n, v, mu, sigma, out_val)
      import :: dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: v
      real(dp), intent(in) :: mu(n)
      real(dp), intent(in) :: sigma
      real(dp), intent(out) :: out_val(n)
    end subroutine pnorm_v

    subroutine plnorm_v(n, v, mu, sigma, sdlog, out_val)
      import :: dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: v
      real(dp), intent(in) :: mu(n)
      real(dp), intent(in) :: sigma
      real(dp), intent(in) :: sdlog
      real(dp), intent(out) :: out_val(n)
    end subroutine plnorm_v

    subroutine pt_v(n, v, mu, sigma, df, out_val)
      import :: dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: v
      real(dp), intent(in) :: mu(n)
      real(dp), intent(in) :: sigma
      real(dp), intent(in) :: df
      real(dp), intent(out) :: out_val(n)
    end subroutine pt_v

    subroutine pgamma_v(n, v, mu, sigma, shape, out_val)
      import :: dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: v
      real(dp), intent(in) :: mu(n)
      real(dp), intent(in) :: sigma
      real(dp), intent(in) :: shape
      real(dp), intent(out) :: out_val(n)
    end subroutine pgamma_v
  end interface

  interface
    ! Sorting function from R
    ! Sorts array x in-place with index tracking
    ! Added: July/2025
    subroutine revsortr(x, idx, n)
      import :: dp
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: idx(n)
      real(dp), intent(inout) :: x(n)
    end subroutine revsortr
  end interface

end module prtree_types
