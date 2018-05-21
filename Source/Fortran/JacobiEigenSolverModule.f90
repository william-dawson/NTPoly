!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing eigenvalues using the Jacobi method.
MODULE JacobiEigenSolverModule
  USE DataTypesModule
  USE DenseMatrixModule
  USE DistributedSparseMatrixModule
  USE DistributedSparseMatrixAlgebraModule
  USE IterativeSolversModule
  USE LoggingModule
  USE MatrixBroadcastModule
  USE MatrixGatherModule
  USE MatrixSendRecvModule
  USE ProcessGridModule
  USE SparseMatrixModule
  USE SparseMatrixAlgebraModule
  USE TripletListModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: JacobiSolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PRIVATE :: SwapData_t
     !> Which process to send the left matrix to.
     INTEGER :: send_left_partner
     !> Which process to send the right matrix to.
     INTEGER :: send_right_partner
     !> A tag to identify where the left matrix is sent to.
     INTEGER :: send_left_tag
     !> A tag to identify where the right matrix is sent to.
     INTEGER :: send_right_tag
     !> Which process to receive the left matrix from.
     INTEGER :: recv_left_partner
     !> Which process to receive the right matrix from.
     INTEGER :: recv_right_partner
     !> A tag to identify the left matrix being received.
     INTEGER :: recv_left_tag
     !> A tag to identify the right matrix being received.
     INTEGER :: recv_right_tag
     !> A full list of the permutation that needs to be performed at each step.
     INTEGER, DIMENSION(:), ALLOCATABLE :: swap_array
  END TYPE SwapData_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PRIVATE :: JacobiData_t
     !! Process Information
     !> Total processors involved in the calculation
     INTEGER :: num_processes
     !> Rank within the 1D process communicator.
     INTEGER :: rank
     !> A communicator for performing the calculation on.
     INTEGER :: communicator
     !! Blocking Information
     INTEGER :: block_start
     INTEGER :: block_end
     INTEGER :: block_rows
     INTEGER :: block_columns
     INTEGER :: rows
     INTEGER :: columns
     INTEGER :: start_row
     INTEGER :: start_column
     INTEGER :: column_divisor
     !! For Inter Process Swaps
     INTEGER, DIMENSION(:), ALLOCATABLE :: phase_array
     !> During the first num_proc rounds
     TYPE(SwapData_t) :: left_swap
     !> After round num_proc+1, we need to do a special swap to maintain the
     !! correct within process order.
     TYPE(SwapData_t) :: mid_swap
     !> In the last num_proc+2->end rounds, we change direction when sending.
     TYPE(SwapData_t) :: right_swap
     !> After the sweeps are finished, this swaps the data back to the original
     !! permutation.
     TYPE(SwapData_t) :: final_swap
  END TYPE JacobiData_t
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using the Block Jacobi Eigenvalue
  !! method.
  !! @param[in] this the matrix to compute the eigendecomposition of.
  !! @param[out] eigenvecotrs the eigenvectors of this.
  !! @param[in] solver_parameters for controlling the solve.
  SUBROUTINE JacobiSolve(this, eigenvectors, solver_parameters)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: solver_parameters
    !! Local Blocking
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: ABlocks
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: VBlocks
    TYPE(SparseMatrix_t) :: local_a
    TYPE(SparseMatrix_t) :: last_a
    TYPE(JacobiData_t) :: jacobi_data
    !! Temporary
    REAL(NTREAL) :: norm_value
    INTEGER :: counter, iteration
    INTEGER :: ierr

    !! Setup Communication
    CALL InitializeJacobi(jacobi_data, this)

    !! Initialize the eigenvectors to the identity.
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, &
         & this%actual_matrix_dimension)
    CALL FillDistributedIdentity(eigenvectors)

    !! Extract to local dense blocks
    ALLOCATE(ABlocks(2,jacobi_data%block_rows))
    CALL GetLocalBlocks(this, jacobi_data, ABlocks)
    ALLOCATE(VBlocks(2,jacobi_data%block_rows))
    CALL GetLocalBlocks(eigenvectors, jacobi_data, VBlocks)
    CALL ComposeSparseMatrix(ABlocks, jacobi_data%block_rows, &
         & jacobi_data%block_columns, local_a)
    DO iteration = 1, solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. iteration .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=iteration-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       !! Loop Over One Jacobi Sweep
       CALL JacobiSweep(ABlocks, VBlocks, jacobi_data, &
            & solver_parameters%threshold)

       !! Compute Norm Value
       CALL CopySparseMatrix(local_a, last_a)
       CALL ComposeSparseMatrix(ABlocks, jacobi_data%block_rows, &
            & jacobi_data%block_columns, local_a)
       CALL IncrementSparseMatrix(local_a,last_a,alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold)
       norm_value = SparseMatrixNorm(last_a)
       CALL MPI_Allreduce(MPI_IN_PLACE, norm_value, 1, MPINTREAL, MPI_SUM,&
            & jacobi_data%communicator, ierr)

       !! Test early exit
       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO

    !! Convert to global matrix
    CALL FillGlobalMatrix(VBlocks, jacobi_data, eigenvectors)

    ! Cleanup
    DO counter = 1, jacobi_data%num_processes*2
       CALL DestructSparseMatrix(ABlocks(1,counter))
       CALL DestructSparseMatrix(ABlocks(2,counter))
       CALL DestructSparseMatrix(VBlocks(1,counter))
       CALL DestructSparseMatrix(VBlocks(2,counter))
    END DO

    DEALLOCATE(ABlocks)
    DEALLOCATE(VBlocks)

    CALL CleanupJacobi(jacobi_data)

  END SUBROUTINE JacobiSolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE InitializeJacobi(jdata, matrix)
    !! Parameters
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: matrix
    !! Local Variables
    INTEGER :: ierr
    !! Music Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: swap0, swap1, swap_temp
    !! Temporary
    INTEGER :: counter
    INTEGER :: stage_counter

    !! Copy The Process Grid Information
    CALL MPI_Comm_dup(matrix%process_grid%within_slice_comm, &
         & jdata%communicator, ierr)
    CALL MPI_Comm_size(jdata%communicator, jdata%num_processes, ierr)
    CALL MPI_Comm_rank(jdata%communicator, jdata%rank, ierr)
    jdata%block_start = (jdata%rank) * 2 + 1
    jdata%block_end = jdata%block_start + 1

    !! Compute the Blocking
    jdata%block_rows = jdata%num_processes*2
    jdata%block_columns = 2

    !! Compute the size of the local matrices
    jdata%rows = matrix%actual_matrix_dimension
    jdata%column_divisor = FLOOR(&
         & 1.0*matrix%actual_matrix_dimension/jdata%num_processes)
    jdata%columns = jdata%column_divisor
    IF (jdata%rank .EQ. jdata%num_processes - 1) THEN
       jdata%columns = matrix%actual_matrix_dimension - &
            & jdata%column_divisor*(jdata%num_processes-1)
    END IF
    jdata%start_column = jdata%column_divisor * jdata%rank + 1
    jdata%start_row = 1

    !! Determine Send Partners.
    ALLOCATE(jdata%phase_array(2*jdata%num_processes-1))
    ALLOCATE(jdata%left_swap%swap_array(2*jdata%num_processes))
    ALLOCATE(jdata%mid_swap%swap_array(2*jdata%num_processes))
    ALLOCATE(jdata%right_swap%swap_array(2*jdata%num_processes))
    ALLOCATE(jdata%final_swap%swap_array(2*jdata%num_processes))
    !! First we create the default permutation.
    ALLOCATE(swap0(2*jdata%num_processes))
    ALLOCATE(swap1(2*jdata%num_processes))
    ALLOCATE(swap_temp(2*jdata%num_processes))
    DO counter = 1, 2*jdata%num_processes
       swap0(counter) = counter
    END DO

    !! Second perform rotations to compute the permutation arrays
    stage_counter = 1
    swap_temp = swap0
    swap1 = swap0
    CALL RotateMusic(jdata,swap1,1)
    jdata%phase_array(stage_counter) = 1
    stage_counter = stage_counter+1
    CALL ComputePartners(jdata, jdata%left_swap, swap_temp, swap1)

    DO counter = 2, jdata%num_processes-1
       CALL RotateMusic(jdata,swap1,1)
       jdata%phase_array(stage_counter) = 1
       stage_counter = stage_counter+1
    END DO

    IF (jdata%num_processes .GT. 1) THEN
       swap_temp = swap1
       CALL RotateMusic(jdata,swap1,2)
       jdata%phase_array(stage_counter) = 2
       stage_counter = stage_counter+1
       CALL ComputePartners(jdata, jdata%mid_swap, swap_temp, swap1)
    END IF

    IF (jdata%num_processes .GT. 2) THEN
       swap_temp = swap1
       CALL RotateMusic(jdata,swap1,3)
       jdata%phase_array(stage_counter) = 3
       stage_counter = stage_counter+1
       CALL ComputePartners(jdata, jdata%right_swap, swap_temp, swap1)
       DO counter = jdata%num_processes+1, 2*jdata%num_processes-3
          CALL RotateMusic(jdata,swap1,3)
          jdata%phase_array(stage_counter) = 3
          stage_counter = stage_counter+1
       END DO
    END IF

    IF (jdata%num_processes .GT. 1) THEN
       CALL ComputePartners(jdata, jdata%final_swap, swap1, swap0)
       jdata%phase_array(stage_counter) = 4
    END IF

    !! Cleanup
    DEALLOCATE(swap0)
    DEALLOCATE(swap1)
    DEALLOCATE(swap_temp)

  END SUBROUTINE InitializeJacobi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct the jacobi data data structure.
  !! @param[inout] jacobi_data the jacobi data.
  SUBROUTINE CleanupJacobi(jacobi_data)
    !! Parameters
    TYPE(JacobiData_t), INTENT(INOUT) :: jacobi_data
    !! Local Variables
    INTEGER :: ierr

    !! MPI Cleanup
    CALL MPI_Comm_free(jacobi_data%communicator, ierr)

    !! Memory Deallocation
    IF (ALLOCATED(jacobi_data%left_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%left_swap%swap_array)
    END IF
    IF (ALLOCATED(jacobi_data%mid_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%mid_swap%swap_array)
    END IF
    IF (ALLOCATED(jacobi_data%right_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%right_swap%swap_array)
    END IF
    IF (ALLOCATED(jacobi_data%final_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%final_swap%swap_array)
    END IF

  END SUBROUTINE CleanupJacobi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetLocalBlocks(distributed, jdata, local)
    !! Parameters
    TYPE(DistributedSparseMatrix_t) :: distributed
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    TYPE(SparseMatrix_t), DIMENSION(:,:) :: local
    !! Local Variables
    TYPE(TripletList_t) :: local_triplets
    TYPE(TripletList_t), DIMENSION(jdata%num_processes) :: send_triplets
    TYPE(TripletList_t) :: received_triplets, sorted_triplets
    !! Some Blocking Information
    INTEGER, DIMENSION(jdata%block_rows) :: divisor
    INTEGER, DIMENSION(2) :: local_divide
    !! Temporary
    TYPE(Triplet_t) :: temp_triplet
    TYPE(SparseMatrix_t) :: local_mat
    INTEGER :: counter, insert
    INTEGER :: ierr

    !! Get The Local Triplets
    CALL ConstructTripletList(local_triplets)
    IF (distributed%process_grid%between_slice_rank .EQ. 0) THEN
       CALL GetTripletList(distributed, local_triplets)
    END IF
    DO counter = 1, jdata%num_processes
       CALL ConstructTripletList(send_triplets(counter))
    END DO
    DO counter = 1, local_triplets%CurrentSize
       CALL GetTripletAt(local_triplets, counter, temp_triplet)
       insert = (temp_triplet%index_column - 1) /  jdata%column_divisor + 1
       IF (insert .GE. jdata%num_processes) THEN
          insert = jdata%num_processes
       END IF
       CALL AppendToTripletList(send_triplets(insert), temp_triplet)
    END DO
    CALL RedistributeTripletLists(send_triplets, &
         & distributed%process_grid%within_slice_comm, received_triplets)
    CALL ShiftTripletList(received_triplets, 0, -(jdata%start_column - 1))
    CALL SortTripletList(received_triplets, jdata%columns, sorted_triplets)
    CALL ConstructFromTripletList(local_mat, sorted_triplets, jdata%rows, &
         & jdata%columns)

    !! Share blocking information
    local_divide(1) = local_mat%columns/2
    local_divide(2) = local_mat%columns - local_divide(1)
    CALL MPI_Allgather(local_divide, 2, MPI_INTEGER, divisor, 2, MPI_INTEGER, &
         & jdata%communicator, ierr)

    !! Split To Blocks
    CALL SplitSparseMatrix(local_mat, jdata%block_rows, jdata%block_columns, &
         & local, block_size_row_in=divisor)

    !! Cleanup
    CALL DestructTripletList(local_triplets)
    DO counter = 1, jdata%num_processes
       CALL DestructTripletList(send_triplets(counter))
    END DO
    CALL DestructTripletList(received_triplets)
    CALL DestructTripletList(sorted_triplets)
    CALL DestructSparseMatrix(local_mat)
  END SUBROUTINE GetLocalBlocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FillGlobalMatrix(local, jdata, global)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(IN) :: local
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: global
    !! Local Variables
    TYPE(SparseMatrix_t) :: TempMat
    TYPE(TripletList_t) :: triplet_list

    !! Get A Global Triplet List and Fill
    CALL ComposeSparseMatrix(local, jdata%block_rows, jdata%block_columns, &
         & TempMat)
    CALL MatrixToTripletList(TempMat, triplet_list)
    CALL ShiftTripletList(triplet_list, jdata%start_row - 1, &
         & jdata%start_column - 1)

    CALL FillFromTripletList(global, triplet_list)

    !! Cleanup
    CALL DestructSparseMatrix(tempmat)
    CALL DestructTripletList(triplet_list)

  END SUBROUTINE FillGlobalMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE JacobiSweep(ABlocks, VBlocks, jdata, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: VBlocks
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(SparseMatrix_t) :: TargetA
    TYPE(SparseMatrix_t) :: TargetV, TargetVT
    INTEGER :: iteration

    !! Loop Over Processors
    DO iteration = 1, jdata%num_processes*2 - 1
       !! Construct A Block To Diagonalize
       CALL ComposeSparseMatrix(ABlocks(:,jdata%block_start:jdata%block_end), &
            & 2, 2, TargetA)

       !! Diagonalize
       CALL DenseEigenDecomposition(TargetA, TargetV, threshold)

       !! Rotation Along Row
       CALL TransposeSparseMatrix(TargetV, TargetVT)
       CALL ApplyToRows(TargetVT, ABlocks, jdata, threshold)

       !! Rotation Along Columns
       CALL ApplyToColumns(TargetV, ABlocks, jdata, threshold)
       CALL ApplyToColumns(TargetV, VBlocks, jdata, threshold)

       !! Swap Blocks
       IF (jdata%num_processes .GT. 1) THEN
          IF (jdata%phase_array(iteration) .EQ. 1) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%left_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%left_swap)
          ELSE IF (jdata%phase_array(iteration) .EQ. 2) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%mid_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%mid_swap)
          ELSE IF (jdata%phase_array(iteration) .EQ. 3) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%right_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%right_swap)
          ELSE
             CALL SwapBlocks(ABlocks, jdata, jdata%final_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%final_swap)
          END IF
       END IF
    END DO

  END SUBROUTINE JacobiSweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Determine the next permutation of rows and columns in a round robin fashion
  !! Based on the music algorithm of \cite{golub2012matrix}.
  !! @param[in] jdata the jacobi data structure.
  !! @param[in] perm_order the current permutation of rows and columns.
  !! @param[in] phase either 1, 2, 3, 4 based on the rotation required.
  PURE SUBROUTINE RotateMusic(jdata, perm_order, phase)
    !! Parameters
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    INTEGER, DIMENSION(2*jdata%num_processes), INTENT(INOUT) :: perm_order
    INTEGER, INTENT(IN) :: phase
    !! Copies of Music
    INTEGER, DIMENSION(jdata%num_processes) :: music_row
    INTEGER, DIMENSION(jdata%num_processes) :: music_column
    INTEGER, DIMENSION(jdata%num_processes) :: music_row_orig
    INTEGER, DIMENSION(jdata%num_processes) :: music_column_orig
    !! Local Variables
    INTEGER :: counter
    INTEGER :: num_pairs
    INTEGER :: ind

    !! For convenience
    num_pairs = jdata%num_processes

    !! First split the permutation order into rows and columns.
    DO counter = 1, num_pairs
       ind = (counter-1)*2 + 1
       music_row(counter) = perm_order(ind)
       music_column(counter) = perm_order(ind+1)
    END DO

    !! Make Copies
    music_row_orig = music_row
    music_column_orig = music_column

    !! No swapping if there is just one processes
    IF (num_pairs .GT. 1) THEN
       IF (phase .EQ. 1) THEN
          !! Rotate Bottom Half
          DO counter = 1, num_pairs - 1
             music_column(counter) = music_column_orig(counter+1)
          END DO
          music_column(num_pairs) = music_row_orig(num_pairs)

          !! Rotate Top Half
          music_row(1) = music_row_orig(1)
          music_row(2) = music_column_orig(1)
          DO counter = 3, num_pairs
             music_row(counter) = music_row_orig(counter-1)
          END DO
       ELSE IF (phase .EQ. 2) THEN
          !! Rotate The Bottom Half
          music_column(1) = music_column_orig(2)
          music_column(2) = music_column_orig(1)
          DO counter = 3, num_pairs
             music_column(counter) = music_row_orig(counter-1)
          END DO

          !! Rotate The Top Half
          music_row(1) = music_row_orig(1)
          DO counter = 2, num_pairs - 1
             music_row(counter) = music_column_orig(counter+1)
          END DO
          music_row(num_pairs) = music_row_orig(num_pairs)
       ELSE IF (phase .EQ. 3) THEN
          !! Rotate Bottom Half
          music_column(1) = music_row_orig(2)
          DO counter = 2, num_pairs
             music_column(counter) = music_column_orig(counter-1)
          END DO

          !! Rotate Top Half
          music_row(1) = music_row_orig(1)
          DO counter = 2, num_pairs - 1
             music_row(counter) = music_row_orig(counter+1)
          END DO
          music_row(num_pairs) = music_column_orig(num_pairs)
       ELSE IF (phase .EQ. 4) THEN
          !! Rotate Bottom Half
          DO counter = 1, num_pairs-1
             music_column(counter) = music_row_orig(counter+1)
          END DO
          music_column(num_pairs) = music_column_orig(num_pairs)

          !! Rotate
          music_row(1) = music_row_orig(1)
          DO counter = 2, num_pairs
             music_row(counter) = music_column_orig(counter-1)
          END DO
       END IF
    END IF

    !! Go from split rows and columns to one big list of pairs.
    DO counter = 1, num_pairs
       ind = (counter-1)*2 + 1
       perm_order(ind) = music_row(counter)
       perm_order(ind+1) = music_column(counter)
    END DO

  END SUBROUTINE RotateMusic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a before and after picture of a permutation, this computes the send
  !! partners required to get this calculation done.
  !! @param[in] jdata the full jacobi_data structure
  !! @param[inout] swap_data a data structure to hold the swap information.
  !! @param[in] perm_before permutation before the swap is done.
  !! @param[in] perm_after permutation after the swap is done.
  PURE SUBROUTINE ComputePartners(jdata, swap_data, perm_before, perm_after)
    !! Parameters
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    TYPE(SwapData_t), INTENT(INOUT) :: swap_data
    INTEGER, DIMENSION(:), INTENT(IN) :: perm_before
    INTEGER, DIMENSION(:), INTENT(IN) :: perm_after
    !! Local Variables
    INTEGER :: send_row
    INTEGER :: send_col
    INTEGER :: recv_row
    INTEGER :: recv_col
    !! Temporary
    INTEGER :: counter
    INTEGER :: inner_counter, outer_counter
    INTEGER :: ind

    !! For simplicitly we extract these into variables
    ind = jdata%rank * 2 + 1
    send_row = perm_before(ind)
    send_col = perm_before(ind+1)
    recv_row = perm_after(ind)
    recv_col = perm_after(ind+1)

    !! Now determine the rank and tag for each of these
    DO counter = 1, 2*jdata%num_processes
       !! Send
       IF (perm_after(counter) .EQ. send_row) THEN
          swap_data%send_left_partner = (counter-1)/jdata%block_columns
          swap_data%send_left_tag = MOD((counter-1), jdata%block_columns)+1
       END IF
       IF (perm_after(counter) .EQ. send_col) THEN
          swap_data%send_right_partner = (counter-1)/jdata%block_columns
          swap_data%send_right_tag = MOD((counter-1), jdata%block_columns)+1
       END IF
       !! Receive
       IF (perm_before(counter) .EQ. recv_row) THEN
          swap_data%recv_left_partner = (counter-1)/jdata%block_columns
       END IF
       IF (perm_before(counter) .EQ. recv_col) THEN
          swap_data%recv_right_partner = (counter-1)/jdata%block_columns
       END IF
    END DO
    swap_data%recv_left_tag = 1
    swap_data%recv_right_tag = 2

    !! Fill In The Swap Data Array For Local Permutation
    DO outer_counter = 1, 2*jdata%num_processes
       DO inner_counter = 1, 2*jdata%num_processes
          IF (perm_after(outer_counter) .EQ. perm_before(inner_counter)) THEN
             swap_data%swap_array(outer_counter) = inner_counter
             EXIT
          END IF
       END DO
    END DO

  END SUBROUTINE ComputePartners
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SwapBlocks(ABlocks, jdata, swap_data)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    TYPE(SwapData_t), INTENT(IN) :: swap_data
    !! Local Matrices
    TYPE(SparseMatrix_t) :: SendLeft, RecvLeft
    TYPE(SparseMatrix_t) :: SendRight, RecvRight
    !! For Sending Data
    TYPE(SendRecvHelper_t) :: send_left_helper
    TYPE(SendRecvHelper_t) :: send_right_helper
    TYPE(SendRecvHelper_t) :: recv_left_helper
    TYPE(SendRecvHelper_t) :: recv_right_helper
    !! Swap Helpers
    TYPE(SparseMatrix_t), DIMENSION(jdata%block_columns,jdata%block_rows) :: &
         & TempABlocks
    INTEGER, DIMENSION(jdata%block_rows) :: split_guide
    !! Stage Tracker
    INTEGER :: completed
    INTEGER :: send_left_stage
    INTEGER :: recv_left_stage
    INTEGER :: send_right_stage
    INTEGER :: recv_right_stage
    INTEGER, PARAMETER :: total_to_complete = 4
    !! Temporary Variables
    INTEGER :: counter
    INTEGER :: index
    INTEGER :: ierr

    !! Swap Rows
    DO counter = 1, jdata%block_rows
       CALL CopySparseMatrix(ABlocks(1,counter), TempABlocks(1,counter))
       CALL CopySparseMatrix(ABlocks(2,counter), TempABlocks(2,counter))
    END DO
    DO counter = 1, jdata%block_rows
       index = swap_data%swap_array(counter)
       CALL CopySparseMatrix(TempABlocks(1,index), ABlocks(1,counter))
       CALL CopySparseMatrix(TempABlocks(2,index), ABlocks(2,counter))
    END DO

    DO counter = 1, jdata%block_rows
       split_guide(counter) = ABlocks(1,counter)%rows
    END DO

    !! Build matrices to swap
    CALL ComposeSparseMatrix(ABlocks(1,:),jdata%block_rows,1,SendLeft)
    CALL ComposeSparseMatrix(ABlocks(2,:),jdata%block_rows,1,SendRight)

    !! Perform Column Swaps
    send_left_stage = 0
    send_right_stage = 0
    recv_left_stage = 0
    recv_right_stage = 0

    completed = 0
    DO WHILE (completed .LT. total_to_complete)
       !! Send Left Matrix
       SELECT CASE(send_left_stage)
       CASE(0) !! Send Sizes
          CALL SendMatrixSizes(SendLeft, swap_data%send_left_partner, &
               & jdata%communicator, send_left_helper, &
               & swap_data%send_left_tag)
          send_left_stage = send_left_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_left_helper)) THEN
             CALL SendMatrixData(SendLeft, swap_data%send_left_partner, &
                  & jdata%communicator, send_left_helper, &
                  & swap_data%send_left_tag)
             send_left_stage = send_left_stage + 1
          END IF
       CASE(2) !! Test Send Outer
          IF (TestSendRecvOuterRequest(send_left_helper)) THEN
             send_left_stage = send_left_stage + 1
          END IF
       CASE(3) !! Test Send Inner
          IF (TestSendRecvInnerRequest(send_left_helper)) THEN
             send_left_stage = send_left_stage + 1
          END IF
       CASE(4) !! Test Send Data
          IF (TestSendRecvDataRequest(send_left_helper)) THEN
             send_left_stage = send_left_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Receive Left Matrix
       SELECT CASE(recv_left_stage)
       CASE(0) !! Receive Sizes
          CALL RecvMatrixSizes(RecvLeft, swap_data%recv_left_partner, &
               & jdata%communicator, recv_left_helper, &
               & swap_data%recv_left_tag)
          recv_left_stage = recv_left_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_left_helper)) THEN
             CALL RecvMatrixData(RecvLeft, swap_data%recv_left_partner, &
                  & jdata%communicator, recv_left_helper, &
                  & swap_data%recv_left_tag)
             recv_left_stage = recv_left_stage + 1
          END IF
       CASE(2) !! Test Receive Outer
          IF (TestSendRecvOuterRequest(recv_left_helper)) THEN
             recv_left_stage = recv_left_stage + 1
          END IF
       CASE(3) !! Test Receive Inner
          IF (TestSendRecvInnerRequest(recv_left_helper)) THEN
             recv_left_stage = recv_left_stage + 1
          END IF
       CASE(4) !! Test Receive Data
          IF (TestSendRecvDataRequest(recv_left_helper)) THEN
             recv_left_stage = recv_left_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Send Right Matrix
       SELECT CASE(send_right_stage)
       CASE(0) !! Send Sizes
          CALL SendMatrixSizes(SendRight, swap_data%send_right_partner, &
               & jdata%communicator, send_right_helper, &
               & swap_data%send_right_tag)
          send_right_stage = send_right_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_right_helper)) THEN
             CALL SendMatrixData(SendRight, swap_data%send_right_partner, &
                  & jdata%communicator, send_right_helper, &
                  & swap_data%send_right_tag)
             send_right_stage = send_right_stage + 1
          END IF
       CASE(2) !! Test Send Outer
          IF (TestSendRecvOuterRequest(send_right_helper)) THEN
             send_right_stage = send_right_stage + 1
          END IF
       CASE(3) !! Test Send Inner
          IF (TestSendRecvInnerRequest(send_right_helper)) THEN
             send_right_stage = send_right_stage + 1
          END IF
       CASE(4) !! Test Send Data
          IF (TestSendRecvDataRequest(send_right_helper)) THEN
             send_right_stage = send_right_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Receive Right Matrix
       SELECT CASE(recv_right_stage)
       CASE(0) !! Receive Sizes
          CALL RecvMatrixSizes(RecvRight, swap_data%recv_right_partner, &
               & jdata%communicator, recv_right_helper, &
               & swap_data%recv_right_tag)
          recv_right_stage = recv_right_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_right_helper)) THEN
             CALL RecvMatrixData(RecvRight, swap_data%recv_right_partner, &
                  & jdata%communicator, recv_right_helper, &
                  & swap_data%recv_right_tag)
             recv_right_stage = recv_right_stage + 1
          END IF
       CASE(2) !! Test Receive Outer
          IF (TestSendRecvOuterRequest(recv_right_helper)) THEN
             recv_right_stage = recv_right_stage + 1
          END IF
       CASE(3) !! Test Receive Inner
          IF (TestSendRecvInnerRequest(recv_right_helper)) THEN
             recv_right_stage = recv_right_stage + 1
          END IF
       CASE(4) !! Test Receive Data
          IF (TestSendRecvDataRequest(recv_right_helper)) THEN
             recv_right_stage = recv_right_stage + 1
             completed = completed + 1
          END IF
       END SELECT
    END DO
    CALL MPI_Barrier(jdata%communicator, ierr)

    CALL SplitSparseMatrix(RecvLeft,jdata%block_rows,1,ABlocks(1,:), &
         & block_size_row_in=split_guide)
    CALL SplitSparseMatrix(RecvRight,jdata%block_rows,1,ABlocks(2,:), &
         & block_size_row_in=split_guide)

    !! Cleanup
    DO counter = 1, jdata%block_rows
       CALL DestructSparseMatrix(TempABlocks(1,counter))
       CALL DestructSparseMatrix(TempABlocks(2,counter))
    END DO
    CALL DestructSparseMatrix(SendLeft)
    CALL DestructSparseMatrix(SendRight)
    CALL DestructSparseMatrix(RecvLeft)
    CALL DestructSparseMatrix(RecvRight)

  END SUBROUTINE SwapBlocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ApplyToColumns(TargetV, ABlocks, jdata, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: TargetV
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Temporary
    TYPE(SparseMatrix_t) :: AMat, TempMat
    INTEGER :: counter, ind
    INTEGER, DIMENSION(2) :: row_sizes, col_sizes

    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       row_sizes(1) = ABlocks(1,ind)%rows
       row_sizes(2) = ABlocks(1,ind+1)%rows
       col_sizes(1) = ABlocks(1,ind)%columns
       col_sizes(2) = ABlocks(2,ind)%columns

       CALL ComposeSparseMatrix(ABlocks(1:2,ind:ind+1), 2, 2, AMat)
       CALL Gemm(AMat, TargetV, TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(:,ind:ind+1), &
            & block_size_row_in=row_sizes, block_size_column_in=col_sizes)
    END DO

    CALL DestructSparseMatrix(AMat)
    CALL DestructSparseMatrix(TempMat)

  END SUBROUTINE ApplyToColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ApplyToRows(TargetV, ABlocks, jdata, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: TargetV
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Temporary
    INTEGER :: counter
    INTEGER :: ind
    TYPE(SparseMatrix_t) :: TempMat
    TYPE(SparseMatrix_t) :: RecvMat
    TYPE(SparseMatrix_t) :: AMat
    INTEGER, DIMENSION(2) :: row_sizes, col_sizes

    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       IF (jdata%rank .EQ. counter-1) THEN
          CALL CopySparseMatrix(TargetV, RecvMat)
       END IF
       CALL BroadcastMatrix(RecvMat, jdata%communicator, counter - 1)

       row_sizes(1) = ABlocks(1,ind)%rows
       row_sizes(2) = ABlocks(1,ind+1)%rows
       col_sizes(1) = ABlocks(1,ind)%columns
       col_sizes(2) = ABlocks(2,ind)%columns

       CALL ComposeSparseMatrix(ABlocks(:,ind:ind+1), 2, 2, AMat)
       CALL Gemm(RecvMat, AMat, TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(:,ind:ind+1), &
            & block_size_row_in=row_sizes, block_size_column_in=col_sizes)
       CALL DestructSparseMatrix(AMat)
    END DO

    CALL DestructSparseMatrix(RecvMat)
    CALL DestructSparseMatrix(TempMat)

  END SUBROUTINE ApplyToRows
END MODULE JacobiEigenSolverModule
