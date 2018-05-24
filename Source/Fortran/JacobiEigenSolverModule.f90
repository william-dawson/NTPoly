!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing eigenvalues using the Jacobi method.
MODULE JacobiEigenSolverModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DenseMatrixModule, ONLY : DenseMatrix_t, ConstructSparseFromDense, &
       & ConstructDenseFromSparse, DenseEigenDecomposition, DestructDenseMatrix
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & ConstructEmptyDistributedSparseMatrix, FillDistributedIdentity, &
       & FillFromTripletList, GetTripletList, CopyDistributedSparseMatrix, &
       & CommSplitDistributedSparseMatrix, DestructDistributedSparseMatrix, &
       & PrintMatrixInformation
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       & PrintIterativeSolverParameters
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader, WriteCitation
  USE MatrixGatherModule, ONLY : BlockingMatrixGather
  USE MatrixSendRecvModule, ONLY : SendRecvHelper_t, RecvSparseMatrixSizes, &
       & RecvSparseMatrixData, SendSparseMatrixSizes, SendSparseMatrixData, &
       & TestSendRecvSizeRequest, TestSendRecvOuterRequest, &
       & TestSendRecvInnerRequest, TestSendRecvDataRequest
  USE SparseMatrixModule, ONLY : SparseMatrix_t, ComposeSparseMatrix, &
       & ConstructFromTripletList, CopySparseMatrix, DestructSparseMatrix, &
       & MatrixToTripletList, SplitSparseMatrix, TransposeSparseMatrix, &
       & PrintSparseMatrix
  USE SparseMatrixAlgebraModule, ONLY : Gemm, IncrementSparseMatrix, &
       & SparseMatrixNorm
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_t, AppendToTripletList, &
       & ConstructTripletList, DestructTripletList, GetTripletAt, &
       & RedistributeTripletLists, ShiftTripletList, SortTripletList
  USE TripletModule, ONLY : Triplet_t
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: JacobiSolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PRIVATE :: SwapData_t
     !> Which process to send the up matrix to.
     INTEGER :: send_up_partner
     !> Which process to send the down matrix to.
     INTEGER :: send_down_partner
     !> A tag to identify where the up matrix is sent to.
     INTEGER :: send_up_tag
     !> A tag to identify where the down matrix is sent to.
     INTEGER :: send_down_tag
     !> Which process to receive the up matrix from.
     INTEGER :: recv_up_partner
     !> Which process to receive the down matrix from.
     INTEGER :: recv_down_partner
     !> A tag to identify the up matrix being received.
     INTEGER :: recv_up_tag
     !> A tag to identify the down matrix being received.
     INTEGER :: recv_down_tag
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
     !> First block row to work on locally..
     INTEGER :: row_block_start
     !> Second block row to work on locally.
     INTEGER :: row_block_end
     !> First block row to work on locally..
     INTEGER :: col_block_start
     !> Second block row to work on locally.
     INTEGER :: col_block_end
     !> How many block rows there are.
     INTEGER :: block_rows
     !> How many block columsn there are.
     INTEGER :: block_columns
     !> Local rows.
     INTEGER :: rows
     !> Local columns.
     INTEGER :: columns
     !> First row stored locally.
     INTEGER :: start_row
     !> First column stored locally.
     INTEGER :: start_column
     !> For evenly dividing the rows.
     INTEGER :: row_divisor
     !! For Inter Process Swaps
     INTEGER, DIMENSION(:), ALLOCATABLE :: phase_array
     !> During the first num_proc rounds
     TYPE(SwapData_t) :: up_swap
     !> After round num_proc+1, we need to do a special swap to maintain the
     !! correct within process order.
     TYPE(SwapData_t) :: mid_swap
     !> In the last num_proc+2->end rounds, we change direction when sending.
     TYPE(SwapData_t) :: down_swap
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
    TYPE(DistributedSparseMatrix_t) :: WorkingMat, WorkingV
    REAL(NTREAL) :: norm_value
    INTEGER :: iteration, II, JJ
    INTEGER :: ierr
    LOGICAL :: active

    !! Write info about the solver
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Eigen Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Jacobi")
       CALL WriteCitation("golub2012matrix")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Handle the case where we have a matrix that is too small
    CALL StartTimer("Fix Size")
    CALL SplitJacobi(this, WorkingMat, active)
    CALL StopTimer("Fix Size")

    IF (active) THEN
       !! Setup Communication
       CALL StartTimer("Initialize-Jacobi")
       CALL InitializeJacobi(jacobi_data, WorkingMat)

       !! Initialize the eigenvectors to the identity.
       CALL ConstructEmptyDistributedSparseMatrix(WorkingV, &
            & this%actual_matrix_dimension, &
            & process_grid_in=WorkingMat%process_grid)
       CALL FillDistributedIdentity(WorkingV)

       !! Extract to local dense blocks
       ALLOCATE(ABlocks(jacobi_data%block_columns,jacobi_data%block_rows))
       CALL GetLocalBlocks(WorkingMat, jacobi_data, ABlocks)
       ALLOCATE(VBlocks(jacobi_data%block_columns,jacobi_data%block_rows))
       CALL GetLocalBlocks(WorkingV, jacobi_data, VBlocks)
       CALL ComposeSparseMatrix(ABlocks, jacobi_data%block_rows, &
            & jacobi_data%block_columns, local_a)
       CALL StopTimer("Initialize-Jacobi")

       !! Main Loop
       CALL StartTimer("Loop-Jacobi")
       DO iteration = 1, solver_parameters%max_iterations
          IF (solver_parameters%be_verbose .AND. iteration .GT. 1) THEN
             CALL WriteListElement(key="Round", int_value_in=iteration-1)
             CALL EnterSubLog
             CALL WriteListElement(key="Convergence", float_value_in=norm_value)
             CALL ExitSubLog
          END IF

          !! Loop Over One Jacobi Sweep
          CALL StartTimer("Sweep-Jacobi")
          CALL JacobiSweep(ABlocks, VBlocks, jacobi_data, &
               & solver_parameters%threshold)
          CALL StopTimer("Sweep-Jacobi")

          !! Compute Norm Value
          CALL CopySparseMatrix(local_a, last_a)
          CALL ComposeSparseMatrix(ABlocks, jacobi_data%block_rows, &
               & jacobi_data%block_columns, local_a)
          CALL IncrementSparseMatrix(local_a,last_a,alpha_in=REAL(-1.0,NTREAL), &
               & threshold_in=solver_parameters%threshold)
          norm_value = SparseMatrixNorm(last_a)
          CALL MPI_Allreduce(MPI_IN_PLACE, norm_value, 1, MPINTREAL, MPI_SUM, &
               & jacobi_data%communicator, ierr)

          !! Test early exit
          IF (norm_value .LE. solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO
       CALL StopTimer("Loop-Jacobi")
    END IF

    !! Convert to global matrix
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, &
         & this%actual_matrix_dimension, &
         & process_grid_in=this%process_grid)
    CALL FillGlobalMatrix(VBlocks, jacobi_data, eigenvectors, active)

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Total_Iterations",int_value_in=iteration)
       CALL PrintMatrixInformation(eigenvectors)
       CALL ExitSubLog
    END IF

    ! Cleanup
    IF (active) THEN
       DO II = 1, jacobi_data%block_rows
          DO JJ = 1, jacobi_data%block_columns
             CALL DestructSparseMatrix(ABlocks(JJ,II))
             CALL DestructSparseMatrix(VBlocks(JJ,II))
          END DO
       END DO

       DEALLOCATE(ABlocks)
       DEALLOCATE(VBlocks)

       CALL DestructDistributedSparseMatrix(WorkingV)
       CALL DestructDistributedSparseMatrix(WorkingMat)

       CALL CleanupJacobi(jacobi_data)
    END IF

  END SUBROUTINE JacobiSolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SplitJacobi(matrix, splitmat, active)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: matrix
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: splitmat
    LOGICAL, INTENT(OUT) :: active
    !! Local Variables
    INTEGER :: mat_dim
    TYPE(DistributedSparseMatrix_t) :: newmat
    INTEGER :: my_color
    LOGICAL :: split_slice

    CALL CopyDistributedSparseMatrix(matrix, splitmat)
    mat_dim = matrix%actual_matrix_dimension
    active = .TRUE.

    DO WHILE(mat_dim/(2*splitmat%process_grid%slice_size) .LT. 1)
       CALL CommSplitDistributedSparseMatrix(splitmat, newmat, my_color, &
            & split_slice)
       CALL CopyDistributedSparseMatrix(newmat, splitmat)
       IF (my_color .GE. 1) THEN
          active = .FALSE.
          EXIT
       END IF
    END DO

    CALL DestructDistributedSparseMatrix(newmat)
  END SUBROUTINE SplitJacobi
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
    jdata%col_block_start = (jdata%rank) * 2 + 1
    jdata%col_block_end = jdata%col_block_start + 1
    jdata%row_block_start = 1
    jdata%row_block_end = 2

    !! Compute the Blocking
    jdata%block_rows = 2
    jdata%block_columns = jdata%num_processes*2

    !! Compute the size of the local matrices
    jdata%columns = matrix%actual_matrix_dimension
    jdata%row_divisor = FLOOR(&
         & 1.0*matrix%actual_matrix_dimension/jdata%num_processes)
    jdata%rows = jdata%row_divisor
    IF (jdata%rank .EQ. jdata%num_processes - 1) THEN
       jdata%rows = matrix%actual_matrix_dimension - &
            & jdata%row_divisor*(jdata%num_processes-1)
    END IF
    jdata%start_row = jdata%row_divisor * jdata%rank + 1
    jdata%start_column = 1

    !! Determine Send Partners.
    ALLOCATE(jdata%phase_array(2*jdata%num_processes-1))
    ALLOCATE(jdata%up_swap%swap_array(2*jdata%num_processes))
    ALLOCATE(jdata%mid_swap%swap_array(2*jdata%num_processes))
    ALLOCATE(jdata%down_swap%swap_array(2*jdata%num_processes))
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
    CALL ComputePartners(jdata, jdata%up_swap, swap_temp, swap1)

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
       CALL ComputePartners(jdata, jdata%down_swap, swap_temp, swap1)
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
    IF (ALLOCATED(jacobi_data%up_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%up_swap%swap_array)
    END IF
    IF (ALLOCATED(jacobi_data%mid_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%mid_swap%swap_array)
    END IF
    IF (ALLOCATED(jacobi_data%down_swap%swap_array)) THEN
       DEALLOCATE(jacobi_data%down_swap%swap_array)
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
    INTEGER, DIMENSION(jdata%block_columns) :: divisor
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
       insert = (temp_triplet%index_row - 1) /  jdata%row_divisor + 1
       IF (insert .GE. jdata%num_processes) THEN
          insert = jdata%num_processes
       END IF
       CALL AppendToTripletList(send_triplets(insert), temp_triplet)
    END DO
    CALL RedistributeTripletLists(send_triplets, &
         & distributed%process_grid%within_slice_comm, received_triplets)
    CALL ShiftTripletList(received_triplets, -(jdata%start_row - 1), 0)
    CALL SortTripletList(received_triplets, jdata%columns, sorted_triplets)
    CALL ConstructFromTripletList(local_mat, sorted_triplets, jdata%rows, &
         & jdata%columns)

    !! Share blocking information
    local_divide(1) = local_mat%rows/2
    local_divide(2) = local_mat%rows - local_divide(1)
    CALL MPI_Allgather(local_divide, 2, MPI_INTEGER, divisor, 2, MPI_INTEGER, &
         & jdata%communicator, ierr)

    !! Split To Blocks
    CALL SplitSparseMatrix(local_mat, jdata%block_rows, jdata%block_columns, &
         & local, block_size_column_in=divisor)

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
  SUBROUTINE FillGlobalMatrix(local, jdata, global, active)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(IN) :: local
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: global
    LOGICAL, INTENT(IN) :: active
    !! Local Variables
    TYPE(SparseMatrix_t) :: TempMat
    TYPE(TripletList_t) :: triplet_list

    !! Get A Global Triplet List and Fill
    IF (active) THEN
       CALL ComposeSparseMatrix(local, jdata%block_rows, jdata%block_columns, &
            & TempMat)
       CALL MatrixToTripletList(TempMat, triplet_list)
       CALL ShiftTripletList(triplet_list, jdata%start_row - 1, &
            & jdata%start_column - 1)
    ELSE
       CALL ConstructTripletList(triplet_list)
    END IF

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
    TYPE(DenseMatrix_t) :: DenseA, DenseV
    INTEGER :: iteration

    !! Loop Over Processors
    DO iteration = 1, jdata%num_processes*2 - 1
       !! Construct A Block To Diagonalize
       CALL ComposeSparseMatrix(ABlocks(&
            & jdata%col_block_start:jdata%col_block_end, &
            & jdata%row_block_start:jdata%row_block_end), &
            & 2, 2, TargetA)

       !! Diagonalize
       CALL StartTimer("Decomposition")
       CALL ConstructDenseFromSparse(TargetA, DenseA)

       CALL DenseEigenDecomposition(DenseA, DenseV)
       CALL ConstructSparseFromDense(DenseV, TargetV)

       CALL DestructDenseMatrix(DenseA)
       CALL DestructDenseMatrix(DenseV)
       CALL StopTimer("Decomposition")

       !! Rotation Along Row
       CALL StartTimer("Row")
       CALL TransposeSparseMatrix(TargetV, TargetVT)
       CALL ApplyToRows(TargetVT, ABlocks, jdata, threshold)
       CALL StopTimer("Row")

       !! Rotation Along Columns
       CALL StartTimer("Column")
       CALL ApplyToColumns(TargetV, ABlocks, jdata, threshold)
       CALL ApplyToColumns(TargetV, VBlocks, jdata, threshold)
       CALL StopTimer("Column")

       !! Swap Blocks
       CALL StartTimer("Swap")
       IF (jdata%num_processes .GT. 1) THEN
          IF (jdata%phase_array(iteration) .EQ. 1) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%up_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%up_swap)
          ELSE IF (jdata%phase_array(iteration) .EQ. 2) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%mid_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%mid_swap)
          ELSE IF (jdata%phase_array(iteration) .EQ. 3) THEN
             CALL SwapBlocks(ABlocks, jdata, jdata%down_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%down_swap)
          ELSE
             CALL SwapBlocks(ABlocks, jdata, jdata%final_swap)
             CALL SwapBlocks(VBlocks, jdata, jdata%final_swap)
          END IF
       END IF
       CALL StopTimer("Swap")
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
          swap_data%send_up_partner = (counter-1)/jdata%block_rows
          swap_data%send_up_tag = MOD((counter-1), jdata%block_rows)+1
       END IF
       IF (perm_after(counter) .EQ. send_col) THEN
          swap_data%send_down_partner = (counter-1)/jdata%block_rows
          swap_data%send_down_tag = MOD((counter-1), jdata%block_rows)+1
       END IF
       !! Receive
       IF (perm_before(counter) .EQ. recv_row) THEN
          swap_data%recv_up_partner = (counter-1)/jdata%block_rows
       END IF
       IF (perm_before(counter) .EQ. recv_col) THEN
          swap_data%recv_down_partner = (counter-1)/jdata%block_rows
       END IF
    END DO
    swap_data%recv_up_tag = 1
    swap_data%recv_down_tag = 2

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
    TYPE(SparseMatrix_t) :: SendUp, RecvUp
    TYPE(SparseMatrix_t) :: SendDown, RecvDown
    !! For Sending Data
    TYPE(SendRecvHelper_t) :: send_up_helper
    TYPE(SendRecvHelper_t) :: send_down_helper
    TYPE(SendRecvHelper_t) :: recv_up_helper
    TYPE(SendRecvHelper_t) :: recv_down_helper
    !! Swap Helpers
    TYPE(SparseMatrix_t), DIMENSION(jdata%block_columns,jdata%block_rows) :: &
         & TempABlocks
    INTEGER, DIMENSION(jdata%block_columns) :: split_guide
    !! Stage Tracker
    INTEGER :: completed
    INTEGER :: send_up_stage
    INTEGER :: recv_up_stage
    INTEGER :: send_down_stage
    INTEGER :: recv_down_stage
    INTEGER, PARAMETER :: total_to_complete = 4
    !! Temporary Variables
    INTEGER :: counter
    INTEGER :: index

    !! Swap Rows
    CALL StartTimer("Swap Rows")
    !$OMP PARALLEL
    !$OMP DO
    DO counter = 1, jdata%block_columns
       CALL CopySparseMatrix(ABlocks(counter,1), TempABlocks(counter,1))
       CALL CopySparseMatrix(ABlocks(counter,2), TempABlocks(counter,2))
    END DO
    !$OMP END DO
    !$OMP DO
    DO counter = 1, jdata%block_columns
       index = swap_data%swap_array(counter)
       CALL CopySparseMatrix(TempABlocks(index,1), ABlocks(counter,1))
       CALL CopySparseMatrix(TempABlocks(index,2), ABlocks(counter,2))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    CALL StopTimer("Swap Rows")

    DO counter = 1, jdata%block_columns
       split_guide(counter) = ABlocks(counter,1)%columns
    END DO

    !! Build matrices to swap
    CALL StartTimer("Compose First")
    CALL ComposeSparseMatrix(ABlocks(:,1),1,jdata%block_columns,SendUp)
    CALL ComposeSparseMatrix(ABlocks(:,2),1,jdata%block_columns,SendDown)
    CALL StopTimer("Compose First")

    !! Perform Column Swaps
    send_up_stage = 0
    send_down_stage = 0
    recv_up_stage = 0
    recv_down_stage = 0

    completed = 0
    CALL StartTimer("Swap Loop")
    DO WHILE (completed .LT. total_to_complete)
       !! Send Up Matrix
       SELECT CASE(send_up_stage)
       CASE(0) !! Send Sizes
          CALL SendSparseMatrixSizes(SendUp, swap_data%send_up_partner, &
               & jdata%communicator, send_up_helper, &
               & swap_data%send_up_tag)
          send_up_stage = send_up_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_up_helper)) THEN
             CALL SendSparseMatrixData(SendUp, swap_data%send_up_partner, &
                  & jdata%communicator, send_up_helper, &
                  & swap_data%send_up_tag)
             send_up_stage = send_up_stage + 1
          END IF
       CASE(2) !! Test Send Outer
          IF (TestSendRecvOuterRequest(send_up_helper)) THEN
             send_up_stage = send_up_stage + 1
          END IF
       CASE(3) !! Test Send Inner
          IF (TestSendRecvInnerRequest(send_up_helper)) THEN
             send_up_stage = send_up_stage + 1
          END IF
       CASE(4) !! Test Send Data
          IF (TestSendRecvDataRequest(send_up_helper)) THEN
             send_up_stage = send_up_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Receive Up Matrix
       SELECT CASE(recv_up_stage)
       CASE(0) !! Receive Sizes
          CALL RecvSparseMatrixSizes(RecvUp, swap_data%recv_up_partner, &
               & jdata%communicator, recv_up_helper, &
               & swap_data%recv_up_tag)
          recv_up_stage = recv_up_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_up_helper)) THEN
             CALL RecvSparseMatrixData(RecvUp, swap_data%recv_up_partner, &
                  & jdata%communicator, recv_up_helper, &
                  & swap_data%recv_up_tag)
             recv_up_stage = recv_up_stage + 1
          END IF
       CASE(2) !! Test Receive Outer
          IF (TestSendRecvOuterRequest(recv_up_helper)) THEN
             recv_up_stage = recv_up_stage + 1
          END IF
       CASE(3) !! Test Receive Inner
          IF (TestSendRecvInnerRequest(recv_up_helper)) THEN
             recv_up_stage = recv_up_stage + 1
          END IF
       CASE(4) !! Test Receive Data
          IF (TestSendRecvDataRequest(recv_up_helper)) THEN
             recv_up_stage = recv_up_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Send Down Matrix
       SELECT CASE(send_down_stage)
       CASE(0) !! Send Sizes
          CALL SendSparseMatrixSizes(SendDown, swap_data%send_down_partner, &
               & jdata%communicator, send_down_helper, &
               & swap_data%send_down_tag)
          send_down_stage = send_down_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_down_helper)) THEN
             CALL SendSparseMatrixData(SendDown, swap_data%send_down_partner, &
                  & jdata%communicator, send_down_helper, &
                  & swap_data%send_down_tag)
             send_down_stage = send_down_stage + 1
          END IF
       CASE(2) !! Test Send Outer
          IF (TestSendRecvOuterRequest(send_down_helper)) THEN
             send_down_stage = send_down_stage + 1
          END IF
       CASE(3) !! Test Send Inner
          IF (TestSendRecvInnerRequest(send_down_helper)) THEN
             send_down_stage = send_down_stage + 1
          END IF
       CASE(4) !! Test Send Data
          IF (TestSendRecvDataRequest(send_down_helper)) THEN
             send_down_stage = send_down_stage + 1
             completed = completed + 1
          END IF
       END SELECT
       !! Receive Down Matrix
       SELECT CASE(recv_down_stage)
       CASE(0) !! Receive Sizes
          CALL RecvSparseMatrixSizes(RecvDown, swap_data%recv_down_partner, &
               & jdata%communicator, recv_down_helper, &
               & swap_data%recv_down_tag)
          recv_down_stage = recv_down_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_down_helper)) THEN
             CALL RecvSparseMatrixData(RecvDown, swap_data%recv_down_partner, &
                  & jdata%communicator, recv_down_helper, &
                  & swap_data%recv_down_tag)
             recv_down_stage = recv_down_stage + 1
          END IF
       CASE(2) !! Test Receive Outer
          IF (TestSendRecvOuterRequest(recv_down_helper)) THEN
             recv_down_stage = recv_down_stage + 1
          END IF
       CASE(3) !! Test Receive Inner
          IF (TestSendRecvInnerRequest(recv_down_helper)) THEN
             recv_down_stage = recv_down_stage + 1
          END IF
       CASE(4) !! Test Receive Data
          IF (TestSendRecvDataRequest(recv_down_helper)) THEN
             recv_down_stage = recv_down_stage + 1
             completed = completed + 1
          END IF
       END SELECT
    END DO
    CALL StopTimer("Swap Loop")

    CALL StartTimer("Split Last")
    CALL SplitSparseMatrix(RecvUp,1,jdata%block_columns,ABlocks(:,1:1), &
         & block_size_column_in=split_guide)
    CALL SplitSparseMatrix(RecvDown,1,jdata%block_columns,ABlocks(:,2:2), &
         & block_size_column_in=split_guide)
    CALL StopTimer("Split Last")

    !! Cleanup
    DO counter = 1, jdata%block_rows
       CALL DestructSparseMatrix(TempABlocks(counter,1))
       CALL DestructSparseMatrix(TempABlocks(counter,2))
    END DO
    CALL DestructSparseMatrix(SendUp)
    CALL DestructSparseMatrix(SendDown)
    CALL DestructSparseMatrix(RecvUp)
    CALL DestructSparseMatrix(RecvDown)

  END SUBROUTINE SwapBlocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ApplyToRows(TargetV, ABlocks, jdata, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: TargetV
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(IN) :: jdata
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Temporary
    TYPE(SparseMatrix_t) :: AMat, TempMat
    INTEGER :: counter, ind
    INTEGER, DIMENSION(2) :: row_sizes, col_sizes

    !$OMP PARALLEL PRIVATE(ind, row_sizes, col_sizes, AMat, TempMat, TempBlocks)
    !$OMP DO
    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       row_sizes(1) = ABlocks(ind,1)%rows
       row_sizes(2) = ABlocks(ind+1,2)%rows
       col_sizes(1) = ABlocks(ind,1)%columns
       col_sizes(2) = ABlocks(ind+1,2)%columns
       CALL ComposeSparseMatrix(ABlocks(ind:ind+1,1:2), 2, 2, AMat)
       CALL Gemm(TargetV, AMat, TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(ind:ind+1,1:2), &
            & block_size_row_in=row_sizes, block_size_column_in=col_sizes)
       CALL DestructSparseMatrix(AMat)
       CALL DestructSparseMatrix(TempMat)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  END SUBROUTINE ApplyToRows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ApplyToColumns(TargetV, ABlocks, jdata, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: TargetV
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    REAL(NTREAL), INTENT(IN) :: threshold
    !! For Receiving
    TYPE(SparseMatrix_t), DIMENSION(jdata%num_processes) :: receive_list
    !! Temporary
    INTEGER :: counter
    INTEGER :: ind
    TYPE(SparseMatrix_t) :: TempMat
    TYPE(SparseMatrix_t) :: AMat
    INTEGER, DIMENSION(2) :: row_sizes, col_sizes

    CALL BlockingMatrixGather(TargetV, receive_list, jdata%communicator)

    !$OMP PARALLEL PRIVATE(ind, row_sizes, col_sizes, AMat, TempMat)
    !$OMP DO
    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       row_sizes(1) = ABlocks(ind,1)%rows
       row_sizes(2) = ABlocks(ind+1,2)%rows
       col_sizes(1) = ABlocks(ind,1)%columns
       col_sizes(2) = ABlocks(ind+1,2)%columns

       CALL ComposeSparseMatrix(ABlocks(ind:ind+1,1:2), 2, 2, AMat)
       CALL Gemm(AMat, receive_list(counter), TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(ind:ind+1,1:2), &
            & block_size_row_in=row_sizes, block_size_column_in=col_sizes)
       CALL DestructSparseMatrix(AMat)
       CALL DestructSparseMatrix(TempMat)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    DO counter = 1, jdata%num_processes
       CALL DestructSparseMatrix(receive_list(counter))
    END DO

    CALL DestructSparseMatrix(TempMat)

  END SUBROUTINE ApplyToColumns
END MODULE JacobiEigenSolverModule
