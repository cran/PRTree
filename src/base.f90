module base
  implicit none

  integer, parameter :: dp = kind(1.d0)

  ! nodes matrix
  type nm
     integer :: var
     real(dp) :: inf
     real(dp) :: sup
     integer :: node
  end type nm

  ! nodes info
  type ni
     integer :: node
     integer :: isTerminal
     integer :: fatherNode
     integer :: depth
     integer :: varCut
     real(dp) :: cutpoint
     logical :: split
  end type ni

  ! tree
  type tr 
     integer :: dim_P(2)
     real(dp), allocatable :: P(:,:)
     real(dp), allocatable :: gamma(:)
     real(dp), allocatable :: yhat(:)
     real(dp) :: MSE
     type(nm), allocatable :: nodes_matrix(:)
     type(ni), allocatable :: nodes_info(:)
     real(dp), allocatable :: sigma(:,:)
     integer :: indep
     integer, allocatable :: XRegion(:)
  end type tr

  interface safe_allocate
     module procedure safe_allocate1
     module procedure safe_allocate2
     module procedure safe_allocate_logical
     module procedure safe_allocate_matrix
     module procedure safe_allocate_matrix_int
     module procedure safe_allocate_real_dp
     module procedure safe_allocate_integer
  end interface safe_allocate

contains

  !--------------------------------
  ! Allocation subroutines
  !--------------------------------
  subroutine safe_allocate1(node, n)
    implicit none
    integer, intent(in) :: n
    type(nm), allocatable, intent(inout) :: node(:)
    if(allocated(node)) deallocate(node)
    allocate(node(n))
    return	
  end subroutine safe_allocate1

  subroutine safe_allocate2(node, n)
    implicit none
    integer, intent(in) :: n
    type(ni), allocatable, intent(inout) :: node(:)
    if(allocated(node)) deallocate(node)
    allocate(node(n))
    return	
  end subroutine safe_allocate2

  subroutine safe_allocate_logical(ind, n)
    implicit none
    integer, intent(in) :: n
    logical, allocatable, intent(inout) :: ind(:)
    if(allocated(ind)) deallocate(ind)
    allocate(ind(n))
    return	
  end subroutine safe_allocate_logical

  subroutine safe_allocate_integer(Z, n)
    implicit none
    integer, intent(in) :: n
    integer, allocatable, intent(inout) :: Z(:)
    if(allocated(Z)) deallocate(Z)
    allocate(Z(n))
    return	
  end subroutine safe_allocate_integer

  subroutine safe_allocate_real_dp(R, n)
    implicit none
    integer, intent(in) :: n
    real(dp), allocatable, intent(inout) :: R(:)
    if(allocated(R)) deallocate(R)
    allocate(R(n))
    return	
  end subroutine safe_allocate_real_dp

  subroutine safe_allocate_matrix(M, nrow, ncol)
    implicit none
    integer, intent(in) :: nrow, ncol
    real(dp), allocatable, intent(inout) :: M(:,:)
    if(allocated(M)) deallocate(M)
    allocate(M(nrow, ncol))
    return	
  end subroutine safe_allocate_matrix

  subroutine safe_allocate_matrix_int(M, nrow, ncol)
    implicit none
    integer, intent(in) :: nrow, ncol
    integer, allocatable, intent(inout) :: M(:,:)
    if(allocated(M)) deallocate(M)
    allocate(M(nrow, ncol))
    return	
  end subroutine safe_allocate_matrix_int

  function inv(n, A) result(Ainv)    
    !-------------------------------------------------------------------------------------------
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    !
    !  DGETRF,  DGETRI
    !-------------------------------------------------------------------------------------------
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n, n)
    real(dp) :: Ainv(n,n)
    real(dp) :: work(n)  ! work array for LAPACK
    integer :: ipiv(n)   ! pivot indices
    integer ::  info

    Ainv = A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    return
  end function inv

  function gammahat(n,k,P,y) result(gamma)
    !---------------------------------------
    !  Calculates gamma = (P'P)^{-1}P'y
    !---------------------------------------
    implicit none
    integer, intent(in) :: n,k
    real(dp), intent(in) :: P(n,k)
    real(dp), intent(in) :: y(n)
    real(dp) :: gamma(k)
    gamma = matmul(matmul(inv(k,matmul(transpose(P), P)),transpose(P)), y)
    return
  end function gammahat

  function pnorm(x, mu, sd) result(fn_val)
    !----------------------------------------------------------
    !  Calculates P(X <= x) for X ~ N(mu, sd^2)
    !----------------------------------------------------------
    implicit none
    real(dp), intent(in) :: x, mu, sd
    real(dp) :: fn_val, pnormr
    fn_val = pnormr(x, mu, sd)
    return
  end function pnorm

  function probs(Iindep, nrow, ncol, X, int_inf, int_sup, sigma, ncol_P) result(prob)
    !--------------------------------------------------------
    ! function to calculate the probabilities in matrix P
    ! using the Gaussian distribution
    !--------------------------------------------------------
    implicit none
    integer, intent(in) :: Iindep
    integer, intent(in) :: nrow, ncol
    real(dp), intent(in) :: X(nrow, ncol)
    real(dp), intent(in) :: int_inf(ncol), int_sup(ncol)
    real(dp), intent(in) :: sigma(max(Iindep*ncol,1), max(Iindep*ncol,1))
    integer, intent(in) :: ncol_P
    real(dp) :: prob(nrow)
    real(dp) :: prob_marg(ncol)
    integer :: i,j

    if(Iindep == 1) then
       ! assuming that the columns in X are independent
       do i = 1,nrow
          ! if any X is missing, use uniform probabilities
          if(any(isnan(X(i,:)))) then
             prob(i) = 1.d0/dble(ncol_P)
          else
             ! use normal distribution with mean X(i,j) 
             do j = 1,ncol
                prob_marg(j) = pnorm(int_sup(j), X(i,j), sigma(1,1)) - &
                     pnorm(int_inf(j), X(i,j), sigma(1,1))
             end do
             prob(i) = product(prob_marg)
          end if
       end do
       return

       ! TO DO:
       ! Multivariate case will be implemented in the future
       
    end if
  end function probs

  subroutine main_calc(y, X, nrow, ncol, Iindep, var_inf, var_sup, sigma,&
       max_terminal_nodes, max_depth, n_min, cp, perc_x, p_min, tree)
    !--------------------------------------------------------------------------
    !  Helper: given y, X and sigma, finds the corresponding tree
    !
    !  Stopping criterias used
    !   - number of terminal nodes (max_terminal_nodes)
    !   - number of observations in the final nodes (n_min)
    !   - node depth (max_depth)
    !   - percentage of probabilities (perc_x) higher than a threshold (p_min)
    !   - reduction in MSE (cp)
    !--------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: Iindep, nrow, ncol
    integer, intent(in) :: max_terminal_nodes, max_depth, n_min
    real(dp), intent(in) :: X(nrow, ncol), y(nrow)
    real(dp), intent(in) :: var_inf(ncol), var_sup(ncol)
    real(dp), intent(in) :: sigma(max((1-Iindep)*ncol,1), max((1-Iindep)*ncol,1))
    real(dp), intent(in) :: cp, perc_x, p_min
    ! The results
    type(tr), intent(out) :: tree

    ! Matrixes of custom types
    type(nm), allocatable :: nodes_matrix(:), temp_nodes_matrix(:)
    type(nm) :: new_nodes(ncol*2), nodes_at(ncol*2)	
    type(ni), allocatable :: nodes_info(:), temp_nodes_info(:)
    type(ni) :: new_nodes_info(2), nodes_info_at(ncol*2) 

    !All the non allocatables
    real(dp) :: MSE_list(2), val1, comp
    real(dp) :: yhat(nrow), MSE, val_at
    integer :: n_terminal_nodes, n_nodes, i, j, k , n_split, node, node_split
    integer :: node_rm, father, var, var_at, ind_Z(1)
    integer :: Xregion(nrow), c1, c2, skip_node, skip_var

    ! All the allocatables
    real(dp), allocatable, dimension(:,:) :: Pant, Pnew
    real(dp), allocatable :: perc_comp(:), gammanew(:)
    integer, allocatable :: terminal_nodes(:), cuttable_nodes(:)
    integer, allocatable :: obs_update(:), obs_update2(:)
    integer, allocatable :: obs_update_at(:), obs_update_at2(:)
    logical, allocatable :: ind(:)

    MSE_list = Huge(1.d0)
    comp = Huge(1.d0)

    !-------------------------------------------------------------------
    ! - The algorithm starts with all data in node 1 (Root)
    ! - The terminal nodes will increase by one in each loop
    ! - The nodes will increase by two in each loop
    ! - The terminal nodes list will be updated at the end of the loop
    !-------------------------------------------------------------------

    !---------------------
    !     Root node
    !---------------------
    Xregion = 1  ! starts with all data in the root

    ! number of terminal nodes and nodes 
    n_terminal_nodes = 1
    n_nodes = 1

    ! position of terminal nodes
    call safe_allocate(terminal_nodes, 1)
    terminal_nodes = 1

    ! Node matrix. In the root node, -inf < X_j < inf
    call safe_allocate(nodes_matrix, ncol)
    nodes_matrix(:)%var = (/(j, j = 1,ncol)/) 
    nodes_matrix(:)%inf = -huge(1.d0)  
    nodes_matrix(:)%sup = Huge(1.d0)          
    nodes_matrix(:)%node = 1

    ! Node info
    call safe_allocate(nodes_info, 1)
    nodes_info(1)%node = 1        
    nodes_info(1)%isTerminal = 1  
    nodes_info(1)%fatherNode = 0  
    nodes_info(1)%depth = 0
    nodes_info(1)%varCut = 0      
    nodes_info(1)%cutpoint = 0.d0
    nodes_info(1)%split = .true.

    !---------------------------------------------------------------------
    ! Loop to create divisions.
    ! Uses a stopping criteria based on the number of terminal nodes
    !---------------------------------------------------------------------
    do while (n_terminal_nodes < max_terminal_nodes)

       !-----------------------------------------------------------------
       ! Stopping criteria:
       !   - checking if there exists at least one node to split
       !-----------------------------------------------------------------
       if(all(nodes_info(terminal_nodes)%split .eqv. .false.)) goto 999

       k = n_terminal_nodes          
       !---------------------------------------------
       ! filling the matrix P with the values 
       !    Psi(x, Rj, sigma), 1 <= j <= ncol,
       ! where Rj is the j-th terminal node (region)
       !---------------------------------------------
       call safe_allocate(Pant, nrow, k) 
       do j = 1,k
          ! find the position of the terminal nodes
          call safe_allocate(ind, size(nodes_matrix%node))
          ind = nodes_matrix%node == terminal_nodes(j)
          ! calculate Psi for each row of X
          Pant(:,j) = probs(Iindep, nrow, ncol, X, pack(nodes_matrix%inf, ind),& 
               pack(nodes_matrix%sup, ind), sigma, k)
       end do

       !-----------------------------------------------------------------
       ! Stopping criteria:
       !  - checking the depth of the nodes and the percentage of
       !    probabilities higher then a given threshold
       !  - only the nodes that pass the test and have split = TRUE
       !    are candidates to split 
       !-----------------------------------------------------------------
       call safe_allocate(perc_comp, k)
       perc_comp =  count(Pant > p_min, dim = 1)/dble(nrow) 

       call safe_allocate(ind, k)
       ind = nodes_info(terminal_nodes)%depth < max_depth .and. perc_comp > perc_x
       ind = ind .and. nodes_info(terminal_nodes)%split  
       if(all(ind .eqv. .false.)) goto 999 

       ! Find the position of the terminal nodes that pass the above criteria
       n_split = count(ind) 
       call safe_allocate(cuttable_nodes, n_split)
       cuttable_nodes = pack(terminal_nodes, ind)               

       !-----------------------------------------------------------------
       ! loop to choose which node to split and which variable to use
       !-----------------------------------------------------------------
       skip_node = 0  
       do node_split = 1, n_split 
          node = cuttable_nodes(node_split)
          skip_var = 0

          do var = 1,ncol   
             ! find the middle point of the interval for the current
             ! variable and node
             ind_Z = findloc(nodes_matrix%var == var .and. &
                  nodes_matrix%node == node, .true.)
             val1 = 0.5*(max(var_inf(var), nodes_matrix(ind_Z(1))%inf) + &
                  min(var_sup(var), nodes_matrix(ind_Z(1))%sup))

             ! skipping criterion:  
             !  - check if the number of observations in the new regions
             !    is too small 
             c1 = count(Xregion == node .and. X(:,var) < val1)
             c2 = count(Xregion == node .and. X(:,var) .ge. val1)

             if(min(c1,c2) < n_min) then
                skip_var = skip_var + 1
                cycle  ! move to the next variable
             end if
             ! if the current node and variable are selected, then the
             ! information on Xregion will be updated for the indexes in obs_update  (NEW)
             call safe_allocate(obs_update, c1)
             obs_update = pack((/(j,j=1,nrow)/), Xregion == node .and. X(:,var) < val1)
             call safe_allocate(obs_update2, c2)
             obs_update2 = pack((/(j,j=1,nrow)/), Xregion == node .and. X(:,var) .ge. val1)

             ! splitting the node using the middle point.
             new_nodes(1:ncol)%node = n_nodes + 1
             new_nodes((ncol+1):2*ncol)%node = n_nodes + 2

             call safe_allocate(ind, size(nodes_matrix%node))
             ind = nodes_matrix%node == node

             new_nodes%sup = (/ pack(nodes_matrix%sup, ind), &
                  pack(nodes_matrix%sup, ind)/)
             new_nodes%inf = (/ pack(nodes_matrix%inf, ind), &
                  pack(nodes_matrix%inf, ind)/)
             new_nodes%var = (/(/(j, j =1,ncol)/), (/(j, j =1,ncol)/)/)
             new_nodes(var)%sup = val1
             new_nodes(var+ncol)%inf = val1

             new_nodes_info%node = (/n_nodes+1, n_nodes+2/)
             new_nodes_info%isTerminal = 1
             new_nodes_info%fatherNode = node
             new_nodes_info%depth = nodes_info(node)%depth+1
             new_nodes_info%varCut = var
             new_nodes_info%cutpoint = val1

             ! creating the new P matrix. Only the last two columns need to be updated
             call safe_allocate(Pnew, nrow, k+1)
             Pnew(:,1:(k-1)) = Pant(:, pack((/(j,j =1,k)/),(/(j,j=1,k)/).ne.node))            
             do j = 1,2            
                ! find the position of the new terminal nodes
                call safe_allocate(ind, size(new_nodes%node))
                ind = new_nodes%node == n_nodes+j
                ! update P 
                Pnew(:,k-1+j) = probs(Iindep, nrow, ncol, X, &
                     pack(new_nodes%inf, ind), pack(new_nodes%sup, ind), sigma, k+1)          
             end do

             ! calculate gamma = (P'P)^{-1}P'y
             call safe_allocate(gammanew, k+1)
             gammanew = gammahat(nrow, k+1, Pnew, y)

             ! checking if this division improves the MSE
             ! If MSE is lower than the previous MSE value then the old values are updated
             MSE = sum((y - MATMUL(Pnew, gammanew))**2)
             if (MSE < comp) then
                ! saving the best split 
                nodes_at = new_nodes
                nodes_info_at = new_nodes_info
                father = node
                var_at = var
                val_at = val1
                comp = MSE
                call safe_allocate(obs_update_at, c1) 
                obs_update_at = obs_update  
                call safe_allocate(obs_update_at2, c2)
                obs_update_at2 = obs_update2  
             end if
          end do

          ! check if any split happened at the current node  
          if(skip_var == ncol) then
             ! from now on, this node wil be skipped since the node cannot
             ! be splitted because at least one region does not have enough data
             nodes_info(node_split)%split = .false.
             skip_node = skip_node + 1
          end if
       end do

       ! check if any split happened  
       if(skip_node == n_split) goto 999

       !-------------------------------------------------
       ! If the code gets here we found a node to split
       ! and the MSE in the node is lower or equal then 
       ! the previous MSE
       !-------------------------------------------------

       ! stopping criteria:
       !  - checking the reduction in MSE. 
       !    Ignore the new node division if the reduction
       !    in the MSE is smaller than cp.
       MSE_list(2) = comp
       if (n_terminal_nodes > 1) then       
          if(1 - MSE_list(2)/MSE_list(1) < cp) goto 999			
       end if
       MSE_list(1) = MSE_list(2)

       !--------------------------------------------------
       ! Updating variables using the new split
       !--------------------------------------------------
       ! temporary variables to save old values
       call safe_allocate(temp_nodes_matrix, ncol*n_nodes)
       call safe_allocate(temp_nodes_info, n_nodes)
       temp_nodes_matrix = nodes_matrix
       temp_nodes_info = nodes_info

       ! Reallocating and saving the updated values
       call safe_allocate(nodes_matrix, ncol*(n_nodes+2))
       call safe_allocate(nodes_info, n_nodes+2)

       ! updating the nodes_matrix
       nodes_matrix(1:ncol*n_nodes) = temp_nodes_matrix
       nodes_matrix((ncol*n_nodes+1):ncol*(n_nodes+2)) = nodes_at

       ! updating the nodes_info
       nodes_info(1:n_nodes) = temp_nodes_info
       nodes_info(n_nodes+1:n_nodes+2) = nodes_info_at	
       nodes_info(father)%isTerminal = 0
       nodes_info(father)%varCut = var_at
       nodes_info(father)%cutpoint = val_at
       nodes_info((n_nodes+1):(n_nodes+2))%split = .true.

       ! Updating the Region information 
       Xregion(obs_update_at) = n_nodes + 1 
       Xregion(obs_update_at2) = n_nodes + 2 

       ! Updating the number of nodes and terminal nodes
       n_terminal_nodes = n_terminal_nodes+1
       n_nodes = n_nodes+2

       ! Updating the list of terminal nodes
       call safe_allocate(ind, n_terminal_nodes)
       ind = nodes_info%isTerminal == 1  
       terminal_nodes = pack(nodes_info%node, ind)	

    end do

999 continue 

    call fill_tree(tree, nrow, ncol, y, X, Iindep, n_terminal_nodes, &
         terminal_nodes, ncol*n_nodes, nodes_matrix, n_nodes, nodes_info, &
         shape(sigma), sigma, Xregion)
    return
  end subroutine main_calc

  subroutine fill_tree(tree, nrow, ncol, y, X, Iindep, n_tn, terminal_nodes, &
       n_no, nodes_matrix, n_inf, nodes_info, dim_s, sigma, Xregion)
    !--------------------------------------------------------------------
    ! Helper: fills the tree variable with the final tree information
    !--------------------------------------------------------------------
    implicit none
    type(tr), intent(inout) :: tree		
    integer, intent(in) :: nrow, ncol, Iindep, n_tn
    integer, intent(in) :: dim_s(2), n_no, n_inf
    integer, intent(in) :: terminal_nodes(n_tn)
    integer, intent(in) :: Xregion(nrow)
    type(nm) :: nodes_matrix(n_no)
    type(ni) :: nodes_info(n_inf)
    real(dp), intent(in) :: y(nrow), X(nrow, ncol)
    real(dp), intent(in) :: sigma(dim_s(1), dim_s(2))		
    real(dp) :: gamma(n_tn), yhat(nrow)
    integer :: j
    logical, allocatable :: ind(:)


    ! The final matrix P
    call safe_allocate(tree%P, nrow, n_tn)
    tree%dim_P = (/nrow, n_tn/)

    !---------------------------------------------
    ! filling the matrix P with the values 
    !    Psi(x, Rj, sigma), 1 <= j <= ncol,
    ! where Rj is the j-th terminal node (region)
    !---------------------------------------------
    do j = 1, n_tn
       ! find the position of the terminal nodes
       call safe_allocate(ind, n_tn)
       ind = nodes_matrix%node == terminal_nodes(j)
       ! calculate Psi for each row of X
       tree%P(:,j) = probs(Iindep, nrow, ncol, X, pack(nodes_matrix%inf, ind), &
            pack(nodes_matrix%sup, ind), sigma, n_tn)
    end do

    ! The final gamma vector and predicted values
    gamma = gammahat(tree%dim_P(1), tree%dim_P(2), tree%P, y)
    yhat =  matmul(tree%P, gamma)
    call safe_allocate(tree%gamma, tree%dim_P(2))
    call safe_allocate(tree%yhat, tree%dim_P(1))
    tree%gamma = gamma	
    tree%yhat = yhat

    ! MSE value for the final tree
    tree%MSE = sum((y-yhat)**2)/dble(tree%dim_P(1))

    ! Final sigma
    call safe_allocate(tree%sigma, dim_s(1), dim_s(2))
    tree%sigma = sigma

    ! Final nodes and info
    call safe_allocate(tree%nodes_matrix, n_no)
    call safe_allocate(tree%nodes_info, n_inf)
    tree%nodes_matrix = nodes_matrix		
    tree%nodes_info = nodes_info

    ! XRegion
    call safe_allocate(tree%XRegion, nrow)
    tree%XRegion = XRegion
    return
  end subroutine fill_tree

  subroutine update_tree(tree, tree_in)
    !--------------------------------------------------------------------
    ! Helper: copy the values from one tree (tree_in) variable
    !         to another (tree)
    !--------------------------------------------------------------------
    implicit none
    type(tr), intent(in) :: tree_in
    type(tr), intent(inout) :: tree
    integer :: dim_nm, dim_ni, dim_sigma(2)

    dim_nm = size(tree_in%nodes_matrix)
    dim_ni = size(tree_in%nodes_info)
    dim_sigma = shape(tree_in%sigma)

    call safe_allocate(tree%P, tree_in%dim_P(1), tree_in%dim_P(2))
    call safe_allocate(tree%gamma, tree_in%dim_P(2))
    call safe_allocate(tree%yhat, tree_in%dim_P(1))
    call safe_allocate(tree%nodes_matrix, dim_nm)
    call safe_allocate(tree%nodes_info, dim_ni)
    call safe_allocate(tree%sigma, dim_sigma(1), dim_sigma(2))
    call safe_allocate(tree%Xregion, tree_in%dim_P(1))

    tree%dim_P = tree_in%dim_P
    tree%P = tree_in%P
    tree%gamma = tree_in%gamma
    tree%yhat = tree_in%yhat 
    tree%MSE = tree_in%MSE 
    tree%nodes_matrix = tree_in%nodes_matrix
    tree%nodes_info = tree_in%nodes_info 
    tree%sigma = tree_in%sigma
    tree%indep = tree_in%indep
    tree%XRegion = tree_in%XRegion
    return
  end subroutine update_tree

  subroutine return_tree(tree, P, dim_P, gamma, yhat, MSE, max_terminal_nodes,&
       nodes_matrix_info, cutpoints, ncol, inf, sup, Xregion)
    !--------------------------------------------------------------
    ! Helper: copy values from final tree (tree) to vectors and
    !         matrices in order to return them to R
    !--------------------------------------------------------------
    implicit none
    type(tr), intent(in) :: tree
    integer, intent(in) :: max_terminal_nodes, ncol
    integer, intent(out) :: nodes_matrix_info(-1+2*max_terminal_nodes, 5)
    integer, intent(inout) :: dim_P(2) 
    integer, intent(inout):: XRegion(tree%dim_P(1))
    real(dp), intent(inout) :: cutpoints(-1+2*max_terminal_nodes)    
    real(dp), intent(inout) :: P(dim_P(1), dim_P(2))
    real(dp), intent(inout) :: gamma(dim_P(2)), yhat(tree%dim_P(1))
    real(dp), intent(inout) :: MSE
    real(dp), intent(inout), dimension(ncol*(-1+2*max_terminal_nodes)) :: inf, sup

    integer :: n_nodes

    ! number of nodes in the final tree
    n_nodes =  size(tree%nodes_info%node)

    ! Updating P 
    dim_P = tree%dim_P
    P(1:dim_P(1), 1:dim_P(2)) = tree%P

    ! Updating gamma
    gamma(1:dim_P(2)) = tree%gamma

    ! predicted values and MSE
    yhat = tree%yhat 
    MSE = tree%MSE

    ! initializing and updating the nodes_matrix_info
    nodes_matrix_info = 0
    nodes_matrix_info(1:n_nodes,1) = tree%nodes_info%node
    nodes_matrix_info(1:n_nodes,2) = tree%nodes_info%isTerminal
    nodes_matrix_info(1:n_nodes,3) = tree%nodes_info%fatherNode
    nodes_matrix_info(1:n_nodes,4) = tree%nodes_info%depth
    nodes_matrix_info(1:n_nodes,5) = tree%nodes_info%varCut

    ! initializing and updating the other variables
    cutpoints = 0.d0
    cutpoints(1:size(tree%nodes_info%cutpoint)) = tree%nodes_info%cutpoint
    inf = 0.d0
    inf(1:ncol*n_nodes) = tree%nodes_matrix%inf
    sup = 0.d0
    sup(1:ncol*n_nodes) = tree%nodes_matrix%sup

    ! updating Xregion
    XRegion = tree%XRegion

    return
  end subroutine return_tree
end module base

subroutine pr_treer(y, X, nrow, ncol, sigma, dim_sigma, &
     max_terminal_nodes, cp, max_depth, n_min, perc_x, p_min, Iindep,&
     P, dim_P, gamma, yhat, MSE, nodes_matrix_info, cutpoints, &
     inf, sup, sigma_best, XRegion)
  use base
  implicit none
  ! R variables
  integer :: nrow, ncol, dim_sigma(2), dim_P(2)
  integer :: max_terminal_nodes, max_depth, Iindep, n_min
  integer :: nodes_matrix_info(-1+2*max_terminal_nodes,5)
  integer :: XRegion(nrow)
  real(dp) :: X(nrow, ncol)
  real(dp) :: y(nrow), yhat(nrow)
  real(dp) :: sigma(dim_sigma(1), dim_sigma(2))  
  real(dp) :: cp, perc_x, p_min  
  real(dp) :: P(dim_P(1), dim_P(2)), gamma(dim_P(2))  
  real(dp) :: MSE
  real(dp) :: cutpoints(-1+2*max_terminal_nodes)
  real(dp), dimension(ncol*(-1+2*max_terminal_nodes)) :: inf, sup
  real(dp) :: sigma_best

  ! auxiliar variables
  real(dp) :: val_at, MSE_temp, gen_dp
  real(dp) :: var_inf(ncol), var_sup(ncol)
  integer :: i, gen_int
  type(tr) :: tree, tree_new

  ! Max and min of each variable 
  var_inf = minval(X, dim = 1, mask = isnan(X) .eqv. .false.)
  var_sup = maxval(X, dim = 1, mask = isnan(X) .eqv. .false.)

  ! initialization
  P = 0.0d0
  gamma = 0.0d0
  sigma_best = 0.d0

  if(Iindep == 0) goto 10	

  MSE_temp = huge(1.d0)
  ! loop to select the best sigma
  do i = 1, dim_sigma(1)	

     call main_calc(y, X, nrow, ncol, Iindep, var_inf, var_sup, sigma(i:i,1:1),  &
          max_terminal_nodes, max_depth, n_min, cp, perc_x, p_min, tree_new) 

     ! For each sigma, the best trees are compared
     ! If necessary, update to the best tree
     if(tree_new%MSE < MSE_temp) then
        MSE_temp = tree_new%MSE
        call update_tree(tree, tree_new)
     end if

  end do

  ! Returning the best sigma
  sigma_best = tree%sigma(1,1)

  goto 20

10 continue	

  ! Find the best tree
  call main_calc(y, X, nrow, ncol, Iindep, var_inf, var_sup, sigma, &
       max_terminal_nodes, max_depth, n_min, cp, perc_x, p_min, tree_new)

20 continue

  ! updating values and returning to R
  call return_tree(tree, P, dim_P, gamma, yhat, MSE, max_terminal_nodes, &
       nodes_matrix_info, cutpoints, ncol, inf, sup, XRegion) 

  return
end subroutine pr_treer


subroutine predict_pr_treer(Iindep, nrow, ncol, X_test, inf, sup, n_terminal_nodes, tn, P, gamma, &
     dim_sigma, sigma, yhat_test)
  use base
  implicit none
  integer :: Iindep, nrow, ncol, n_terminal_nodes 
  integer :: tn(n_terminal_nodes), dim_sigma(2)
  real(dp) :: X_test(nrow, ncol), gamma(n_terminal_nodes), yhat_test(nrow)
  real(dp) :: sigma(dim_sigma(1), dim_sigma(2))  
  real(dp), dimension(ncol*(-1+2*n_terminal_nodes)) :: inf, sup
  real(dp) :: P(nrow, n_terminal_nodes)
  real(dp) :: infj(ncol), supj(ncol)
  integer :: j

  !---------------------------------------------
  ! filling the matrix P with the values 
  !    Psi(x, Rj, sigma), 1 <= j <= ncol,
  ! where Rj is the j-th terminal node (region)
  !---------------------------------------------
  do j = 1, n_terminal_nodes            
     ! calculate Psi for each row of X
     infj = inf((ncol*(tn(j)-1)+1):ncol*(tn(j)))
     supj = sup((ncol*(tn(j)-1)+1):ncol*(tn(j)))
     P(:,j) = probs(Iindep, nrow, ncol, X_test, &
          infj, supj, sigma, n_terminal_nodes)       
  end do
  yhat_test = matmul(P, gamma)
  return
end subroutine predict_pr_treer

