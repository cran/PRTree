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
  real(dp), parameter :: pos_inf = huge(1.0_dp)          ! Positive
  real(dp), parameter :: neg_inf = -huge(1.0_dp)         ! Negative
  real(dp), parameter :: eps = epsilon(1.0_dp)           ! Machine epsilon

  ! Threshold for feature importance
  real(dp), parameter :: feature_threshold = 1.0e-7_dp

  ! Maximum depth for direct probability calculation 1/2^depth.
  ! This is derived from machine epsilon to avoid magic numbers. A depth
  ! greater than this results in a probability smaller than epsilon(1.0_dp),
  ! which is numerically negligible. For double precision, this is 53.
  integer, parameter :: max_prob_depth = exponent(1.0_dp / eps) + 1

  ! Logging levels for printInfo
  integer :: printinfo = -1
  integer, parameter :: log_quiet = -1              ! No output
  integer, parameter :: log_base = 0                ! Basic output (errors, final results)
  integer, parameter :: log_verbose = 1             ! More detailed info on splits
  integer, parameter :: log_detailed = 2            ! Verbose info on candidates and MSE
  integer, parameter :: log_debug = 3               ! Full debug output (indexes, thresholds)
  integer, parameter :: log_debug_deep = 10         ! Deep debug

  ! ids for left and righ nodes. Do not change this values        !
  integer, parameter :: left_id = 1          ! id for the left node
  integer, parameter :: right_id = 2         ! id for the right node

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
  !  = > tree_ctrl: tree building criteria
  !  = > tree_net: tree nodes, matrices and vectors
  !    == > tree_node: individual tree nodes
  !        == = > node_state: state information for each node
  !            == > info_split: information on splitting
  !=================================
  type idx_prob
    !---------------------------------------------
    ! helper variable to simplify functions calls
    ! Added July/2025
    !---------------------------------------------
    integer :: n_miss(2) = 0
    integer :: n_prob = 0    ! number of indexes with P > eps
    integer, allocatable :: idx(:)       ! idx completo and P > eps
    integer, allocatable :: idx_na_f(:)  ! idx Xf miss and P > eps
    integer, allocatable :: idx_any_na(:)! idx X miss and P > eps
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
    real(dp), allocatable :: p(:, :) ! Probability of the new nodes
    logical, allocatable :: na_info(:, :) ! mask for NA and P > eps
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

  type node_state
    !---------------------------------------------------------------------------
    ! Node statistics. Information to help splitting the node
    !  - use one for each node
    !  - deallocate after splitting to save memory
    ! Added: July/2025
    !---------------------------------------------------------------------------
    ! constants
    logical :: update = .true.          ! if the node needs update
    integer :: n_cand_found = 0         ! number of candidates found in this node
    integer :: n_obs = 0                ! node size
    integer :: n_imp                    ! number of important features for this node
    !
    ! allocatable arrays: update during the split process
    !
    integer, allocatable :: idx(:)               ! indexes of X in this node
    integer, allocatable :: imp_feat(:)          ! important features for this node
    type(info_thr), allocatable :: thresholds(:) ! thresholds for each variable
    logical, allocatable :: na_info(:, :)        ! na mask for X_f (P > eps) and/or X (P > eps)
    !
    ! splitting candidates: computed only once for each node
    !
    type(info_split), allocatable :: best(:) ! splitting candidates
  end type node_state

  type tree_node
    !----------------------------------------------------------------------
    ! Node information.
    ! Use one for each node in the tree
    ! Modified July/2025
    !----------------------------------------------------------------------
    ! constants
    integer :: id = 0              ! node id
    integer :: isterminal = 1      ! leaf = 1, otherwhise = 0
    integer :: fathernode = 0      ! from which node it originated
    integer :: depth = 0           ! depth in the tree
    integer :: feature = -1        ! variable that defines the split
    real(dp) :: threshold = 0.0_dp ! threshold
    logical :: split = .true.      ! splitting condition
    !
    ! node information: Free (only) state after splitting
    ! bounds: do not deallocate. Will be used to update the bounds matrix
    !
    real(dp), allocatable :: bounds(:, :) ! region bounds (n_feat x 2)
    type(node_state) :: state  ! splitting information
  end type tree_node

  type tree_net
    !------------------------------------------------------------------
    ! Tree nodes and the corresponding parameters and predicitions
    ! Added July/2025
    !------------------------------------------------------------------
    ! constants
    integer :: n_nodes = 1      ! number of nodes
    integer :: n_tn = 1         ! number of terminal nodes
    integer :: n_cand_found = 0 ! number of candidates found
    !
    ! allocatable arrays: update after each split
    !
    integer, allocatable :: region(:)    ! region for each X (n_obs)
    integer, allocatable :: tn_id(:)     ! id for terminal node (n_tn)
    integer, allocatable :: n_zero(:)    ! number of near zero probabilities (n_tn)
    logical, allocatable :: p_zero(:, :) ! mask for P <= eps (n_obs x n_tn)
    real(dp), allocatable :: p(:, :)     ! matrix P (n_obs x n_tn)
    real(dp), allocatable :: gammahat(:)    ! coefficients (n_tn)
    real(dp), allocatable :: yhat(:)     ! predictions (n_obs)
    real(dp) :: mse = pos_inf            ! mean squared errors
    !
    ! nodes in the tree: update after each split
    !
    type(tree_node), allocatable :: nodes(:) ! nodes object (n_nodes)
  end type tree_net

  type tree_ctrl
    !-----------------------------------------------------------------------
    ! criteria used to build the tree
    ! Added July/2025
    !-------------------------------------------------------------------
    integer :: crit = 3                    ! criterion to assign missing
    integer :: max_tn = 1                  ! Max terminal nodes
    integer :: max_d = 1                   ! Max depth
    integer :: min_obs = 1                 ! Min observations
    integer :: n_cand = 3                  ! Candidates count
    logical :: by_node = .false.           ! find best candidate by node or global
    real(dp) :: min_prop = 0.1_dp          ! Min proportion
    real(dp) :: min_prob = 0.05_dp         ! Min probability
    real(dp) :: cp = 0.01_dp               ! Complexity param (reduction in mse)
  end type tree_ctrl

  type tree_data
    !------------------------------------------------------------------
    ! information on the input data
    ! Added July/2025
    !------------------------------------------------------------------
    ! constants
    integer :: n_obs = 1  ! observations
    integer :: n_feat = 1 ! features
    integer :: fill = 2   ! method used to fill P
    !
    ! input data
    !
    real(dp), allocatable :: y(:)     ! Response (n_obs)
    real(dp), allocatable :: x(:, :)  ! covariates (n_obs x n_feat)
    real(dp), allocatable :: sigma(:) ! distribution parameter (n_feat)
    !
    ! allocatable arrays: for missing data manipulation
    !
    logical, allocatable :: na(:, :)  ! missing flag (n_obs x n_feat)
    integer, allocatable :: n_miss(:) ! number of missing cases for each feature
    logical :: any_na = .false. ! if the data has any NA
  end type tree_data

  type tree_tdata
    !------------------------------------------------------------------
    ! information on the input data
    ! Added July/2025
    !------------------------------------------------------------------
    ! constants
    integer :: n_obs = 1          ! observations
    real(dp), allocatable :: x(:, :)         ! covariates (n_obs x n_feat)
    real(dp), allocatable :: y(:)            ! Response (n_obs
    real(dp), allocatable :: yhat(:)         ! predicitons for the response (n_obs)
    real(dp), allocatable :: p(:, :)         ! probability matrix the test data
    real(dp) :: mse = pos_inf                ! mean squared errors for test data
  end type tree_tdata

  type argsdist
    !-------------------------------------------------------------
    ! Distribution related parameters
    ! Added July/2025
    !-------------------------------------------------------------
    integer :: dist_id = 1           ! code identifying the distribution
    real(dp) :: dist_par = 0.0_dp         ! extra parameter (if any)
  end type argsdist

  type tree_model
    !-------------------------------------------------------------
    ! tree object
    ! Added July/2025
    !-------------------------------------------------------------
    type(tree_data) :: train         ! information on the training data
    type(tree_tdata) :: test         ! information on the testing data
    type(tree_ctrl) :: ctrl          ! tree building criteria
    type(tree_net) :: net            ! tree nodes, matrices and vectors
    type(argsdist) :: dist           ! distribution related parameters
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
    integer :: n_start = 0   ! threshold search
    integer :: direction = 1 ! threshold search
  end type h_data

  !=================================
  ! Interfaces to R functions
  !=================================
  interface
    !--------------------------------------------------------------------
    ! Interfaces for external C functions from R's math library.
    ! By declaring these functions as 'pure', we assert to the compiler
    ! that they have no side effects. This is necessary to allow them
    ! to be called from 'elemental' (and therefore 'pure') Fortran
    ! procedures.
    !--------------------------------------------------------------------
    ! Added: July/2025
    pure function pnorm_pure(x, mu, sigma) result(res)
      import :: dp
      implicit none
      real(dp), intent(in) :: x, mu, sigma
      real(dp) :: res
    end function pnorm_pure

    pure function plnorm_pure(x, meanlog, sdlog) result(res)
      import :: dp
      implicit none
      real(dp), intent(in) :: x, meanlog, sdlog
      real(dp) :: res
    end function plnorm_pure

    pure function pt_pure(x, df) result(res)
      import :: dp
      implicit none
      real(dp), intent(in) :: x, df
      real(dp) :: res
    end function pt_pure

    pure function pgamma_pure(x, shape, scale) result(res)
      import :: dp
      implicit none
      real(dp), intent(in) :: x, shape, scale
      real(dp) :: res
    end function pgamma_pure
  end interface

  interface
    ! Sorting function from R
    ! Sorts array x in-place with index tracking
    ! Added: July/2025
    pure subroutine revsortr(x, idx, n)
      import :: dp
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: idx(n)
      real(dp), intent(inout) :: x(n)
    end subroutine revsortr
  end interface

  interface vector
    ! These functions fill the output with the value in the input.
    ! To be used to avoid calling allocate
    ! Added July/2025
    module procedure vector_int
    module procedure vector_real
    module procedure vector_threshold
    module procedure vector_split
  end interface vector

  interface matrix
    ! These functions fill the output with the value in the input.
    ! To be used to avoid calling allocate
    ! Added July/2025
    module procedure matrix_int
    module procedure matrix_real
    module procedure matrix_logical
  end interface matrix

contains
  pure function vector_int(n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates an integer vector of length n filled with source.
    !
    ! ARGUMENTS
    !   n      [in] : Integer, length of the vector.
    !   source [in] : Integer, value to fill the vector.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: source
    integer :: f(n)
    f = source
  end function vector_int

  pure function vector_real(n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a real(dp) vector of length n filled with source.
    !
    ! ARGUMENTS
    !   n      [in] : Integer, length of the vector.
    !   source [in] : Real(dp), value to fill the vector.
    !
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: source
    real(dp) :: f(n)
    f = source
  end function vector_real

  pure function vector_threshold(n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a type(info_thr) vector of length n filled with thr.
    !
    ! ARGUMENTS
    !   n      [in] : Integer, length of the vector.
    !   source [in] : Type(info_thr), value to fill the vector.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    type(info_thr), intent(in) :: source
    type(info_thr) :: f(n)
    f = source
  end function vector_threshold

  pure function vector_split(n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a type(info_split) vector of length n filled with source.
    !
    ! ARGUMENTS
    !   n      [in] : Integer, length of the vector.
    !   source [in] : Type(info_split), value to fill the vector.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    type(info_split), intent(in) :: source
    type(info_split) :: f(n)
    f = source
  end function vector_split

  pure function matrix_int(m, n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates an integer matrix of size m x n filled with x.
    !
    ! ARGUMENTS
    !   m      [in] : Integer, number of rows.
    !   n      [in] : Integer, number of columns.
    !   source [in] : Integer, value to fill the matrix.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, source
    integer :: f(m, n)
    f = source
  end function matrix_int

  pure function matrix_real(m, n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a real(dp) matrix of size m x n filled with x.
    !
    ! ARGUMENTS
    !   m      [in] : Integer, number of rows.
    !   n      [in] : Integer, number of columns.
    !   source [in] : Real(dp), value to fill the matrix.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n
    real(dp), intent(in) :: source
    real(dp) :: f(m, n)
    f = source
  end function matrix_real

  pure function matrix_logical(m, n, source) result(f)
    !----------------------------------------------------------------------------
    ! Creates a logical matrix of size m x n filled with x.
    !
    ! ARGUMENTS
    !   m      [in] : Integer, number of rows.
    !   n      [in] : Integer, number of columns.
    !   source [in] : logical, value to fill the matrix.
    !
    ! Added July/2025
    !----------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n
    logical, intent(in) :: source
    logical :: f(m, n)
    f = source
  end function matrix_logical

end module prtree_types
