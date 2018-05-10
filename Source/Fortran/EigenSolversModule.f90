!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing eigenvalues
MODULE EigenSolversModule
  USE DataTypesModule
  USE DenseMatrixModule
  USE DistributedSparseMatrixModule
  USE DistributedSparseMatrixAlgebraModule
  USE IterativeSolversModule
  USE LoggingModule
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
  PUBLIC :: DistributedEigenDecomposition
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
     INTEGER :: block_dimension
     INTEGER :: block_rows
     INTEGER :: block_columns
     INTEGER :: rows
     INTEGER :: columns
     INTEGER :: start_row
     INTEGER :: start_column
     !! For Inter Process Swaps
     INTEGER :: send_left_partner
     INTEGER :: send_right_partner
     INTEGER :: send_left_tag
     INTEGER :: send_right_tag
     INTEGER :: recv_left_partner
     INTEGER :: recv_right_partner
     INTEGER :: recv_left_tag
     INTEGER :: recv_right_tag
     !! For Keeping Track of Sweeps
     ! INTEGER, DIMENSION(:), ALLOCATABLE :: initial_music_row
     ! INTEGER, DIMENSION(:), ALLOCATABLE :: initial_music_column
     ! INTEGER, DIMENSION(:), ALLOCATABLE :: music_row
     ! INTEGER, DIMENSION(:), ALLOCATABLE :: music_column
     INTEGER, DIMENSION(:), ALLOCATABLE :: music_swap
  END TYPE JacobiData_t
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DistributedEigenDecomposition(this, eigenvectors, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Blocking
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: ABlocks
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: VBlocks
    TYPE(SparseMatrix_t) :: local_v
    TYPE(SparseMatrix_t) :: last_v
    TYPE(JacobiData_t) :: jacobi_data
    !! Temporary
    REAL(NTREAL) :: norm_value
    INTEGER :: counter, iteration

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    !! Setup Communication
    CALL InitializeJacobi(jacobi_data, this)

    !! Initialize the eigenvectors to the identity.
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, &
         & this%actual_matrix_dimension)
    CALL FillDistributedIdentity(eigenvectors)

    !! Extract to local dense blocks
    ALLOCATE(ABlocks(2,slice_size*2))
    CALL GetLocalBlocks(this, jacobi_data, ABlocks)
    ALLOCATE(VBlocks(2,slice_size*2))
    CALL GetLocalBlocks(eigenvectors, jacobi_data, VBlocks)
    CALL ComposeSparseMatrix(VBlocks, jacobi_data%block_rows, &
         & jacobi_data%block_columns, local_v)

    ! DO iteration = 1, solver_parameters%max_iterations
    DO iteration = 1, 2
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
       CALL CopySparseMatrix(local_v, last_v)
       CALL ComposeSparseMatrix(VBlocks, jacobi_data%block_rows, &
            & jacobi_data%block_columns, local_v)
       CALL IncrementSparseMatrix(local_v,last_v,alpha_in=REAL(-1.0,NTREAL), &
            & threshold_in=solver_parameters%threshold)
       norm_value = SparseMatrixNorm(last_v)

       !! Test early exit
       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO

    !! Convert to global matrix
    CALL FillGlobalMatrix(VBlocks, jacobi_data, eigenvectors)

    !! Cleanup
    DO counter = 1, jacobi_data%num_processes*2
       CALL DestructSparseMatrix(ABlocks(1,counter))
       CALL DestructSparseMatrix(ABlocks(2,counter))
       CALL DestructSparseMatrix(VBlocks(1,counter))
       CALL DestructSparseMatrix(VBlocks(2,counter))
    END DO

    DEALLOCATE(ABlocks)
    DEALLOCATE(VBlocks)

    CALL CleanupJacobi(jacobi_data)

  END SUBROUTINE DistributedEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE InitializeJacobi(jdata, matrix)
    !! Parameters
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: matrix
    !! Local Variables
    INTEGER :: matrix_dimension
    INTEGER :: iteration
    INTEGER :: ierr
    !! Music Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: music_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: music_column

    !! Copy The Process Grid Information
    jdata%num_processes = slice_size
    jdata%rank = within_slice_rank
    CALL MPI_Comm_dup(within_slice_comm, jdata%communicator, ierr)
    jdata%block_start = (jdata%rank) * 2 + 1
    jdata%block_end = jdata%block_start + 1

    !! Initialize the Music Arrays
    ALLOCATE(music_row(jdata%num_processes))
    ALLOCATE(music_column(jdata%num_processes))
    DO iteration = 1, slice_size
       music_row(iteration) = (iteration-1)*2 + 1
       music_column(iteration) = (iteration-1)*2 + 2
    END DO

    !! Determine swap partners for music by doing one rotation
    CALL RotateMusic(jdata, music_row, music_column)
    ALLOCATE(jdata%music_swap(2*jdata%num_processes))
    DO iteration = 1, jdata%num_processes
       jdata%music_swap(2*iteration-1) = music_row(iteration)
       jdata%music_swap(2*iteration) = music_column(iteration)
    END DO

    !! Determine Send Partners
    jdata%send_left_partner = jdata%rank + 1
    jdata%send_left_tag = 1
    IF (jdata%rank .EQ. 0) THEN
       jdata%send_left_partner = 0
    ELSE IF (jdata%rank .EQ. jdata%num_processes - 1) THEN
       jdata%send_left_partner = slice_size - 1
       jdata%send_left_tag = 2
    END IF
    jdata%send_right_partner = jdata%rank - 1
    jdata%send_right_tag = 2
    IF (jdata%num_processes .EQ. 1) THEN
       jdata%send_right_partner = 0
       jdata%send_right_tag = 2
    ELSE IF (jdata%rank .EQ. 0) THEN
       jdata%send_right_partner = jdata%rank + 1
       jdata%send_right_tag = 1
    END IF

    !! Determine Recv Partners
    jdata%recv_left_partner = jdata%rank - 1
    jdata%recv_left_tag = 1
    IF (jdata%rank .EQ. 0) THEN
       jdata%recv_left_partner = 0
    END IF
    jdata%recv_right_partner = jdata%rank + 1
    jdata%recv_right_tag = 2
    IF (jdata%rank .EQ. jdata%num_processes - 1) THEN
       jdata%recv_right_partner = jdata%num_processes - 1
    END IF

    !! Compute the Blocking
    matrix_dimension = matrix%actual_matrix_dimension
    jdata%block_rows = jdata%num_processes*2
    jdata%block_columns = 2
    jdata%block_dimension = CEILING(matrix_dimension/(1.0*jdata%block_rows))
    jdata%rows = jdata%block_dimension * jdata%block_rows
    jdata%columns = jdata%block_dimension * jdata%block_columns
    jdata%start_column = jdata%columns * within_slice_rank + 1
    jdata%start_row = 1

    !! Cleanup
    DEALLOCATE(music_row)
    DEALLOCATE(music_column)

  END SUBROUTINE InitializeJacobi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CleanupJacobi(jacobi_data)
    !! Parameters
    TYPE(JacobiData_t), INTENT(INOUT) :: jacobi_data
    !! Local Variables
    INTEGER :: ierr

    CALL MPI_Comm_free(jacobi_data%communicator, ierr)

    IF (ALLOCATED(jacobi_data%music_swap)) THEN
       DEALLOCATE(jacobi_data%music_swap)
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
    TYPE(TripletList_t), DIMENSION(slice_size) :: send_triplets
    TYPE(TripletList_t) :: received_triplets, sorted_triplets
    !! Temporary
    TYPE(Triplet_t) :: temp_triplet
    TYPE(SparseMatrix_t) :: local_mat
    INTEGER :: counter, insert

    !! Get The Local Triplets
    CALL GetTripletList(distributed, local_triplets)
    DO counter = 1, jdata%num_processes
       CALL ConstructTripletList(send_triplets(counter))
    END DO
    DO counter = 1, local_triplets%CurrentSize
       CALL GetTripletAt(local_triplets, counter, temp_triplet)
       insert = (temp_triplet%index_column - 1) / jdata%columns + 1
       CALL AppendToTripletList(send_triplets(insert), temp_triplet)
    END DO
    CALL RedistributeTripletLists(send_triplets, within_slice_comm, &
         & received_triplets)
    CALL ShiftTripletList(received_triplets, 0, -(jdata%start_column - 1))
    CALL SortTripletList(received_triplets, jdata%columns, sorted_triplets)
    CALL ConstructFromTripletList(local_mat, sorted_triplets, jdata%rows, &
         & jdata%columns)

    !! Split To Blocks
    CALL SplitSparseMatrix(local_mat, jdata%block_rows, jdata%block_columns, &
         & local)

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
       CALL SwapBlocks(ABlocks, jdata)
       
    END DO

    ! CALL MPI_Barrier(global_comm, grid_error)
    ! IF (global_rank .EQ. 1) THEN
    !    ! WRITE(*,*) "RANK 1"
    !    CALL PrintSparseMatrix(ABlocks(1,jdata%block_start))
    !    CALL PrintSparseMatrix(ABlocks(2,jdata%block_start))
    !    CALL PrintSparseMatrix(ABlocks(1,jdata%block_end))
    !    CALL PrintSparseMatrix(ABlocks(2,jdata%block_end))
    !    WRITE(*,*)
    ! END IF
    ! CALL MPI_Barrier(global_comm, grid_error)

  END SUBROUTINE JacobiSweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RotateMusic(jdata, music_row, music_column)
    !! Parameters
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    INTEGER, DIMENSION(jdata%num_processes) :: music_row
    INTEGER, DIMENSION(jdata%num_processes) :: music_column
    !! Copies of Music
    INTEGER, DIMENSION(jdata%num_processes) :: music_row_orig
    INTEGER, DIMENSION(jdata%num_processes) :: music_column_orig
    !! Local Variables
    INTEGER :: counter
    INTEGER :: num_pairs

    !! No swapping if there is just one processes
    IF (jdata%num_processes .GT. 1) THEN
       !! For convenience
       num_pairs = jdata%num_processes

       !! Make Copies
       music_row_orig = music_row
       music_column_orig = music_column

       !! Rotate Bottom Half
       music_column(num_pairs) = music_row_orig(num_pairs)
       DO counter = 1, num_pairs - 1
          music_column(counter) = music_column_orig(counter+1)
       END DO

       !! Rotate Top Half
       IF(num_pairs .GT. 1) THEN
          music_row(2) = music_column_orig(1)
       END IF
       DO counter = 3, num_pairs
          music_row(counter) = music_row_orig(counter-1)
       END DO
    END IF

  END SUBROUTINE RotateMusic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SwapBlocks(ABlocks, jdata)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:,:), INTENT(INOUT) :: ABlocks
    TYPE(JacobiData_t), INTENT(INOUT) :: jdata
    !! Local Matrices
    TYPE(SparseMatrix_t) :: SendLeft, RecvLeft
    TYPE(SparseMatrix_t) :: SendRight, RecvRight
    TYPE(SparseMatrix_t), DIMENSION(jdata%block_columns,jdata%block_rows) :: &
         & TempABlocks
    !! For Sending Data
    TYPE(SendRecvHelper_t) :: send_left_helper
    TYPE(SendRecvHelper_t) :: send_right_helper
    TYPE(SendRecvHelper_t) :: recv_left_helper
    TYPE(SendRecvHelper_t) :: recv_right_helper
    !! Temporary Variables
    INTEGER :: completed
    INTEGER :: send_left_stage
    INTEGER :: recv_left_stage
    INTEGER :: send_right_stage
    INTEGER :: recv_right_stage
    INTEGER, PARAMETER :: total_to_complete = 4
    INTEGER :: counter
    INTEGER :: index
    INTEGER :: ierr

    !! Swap Rows
    DO counter = 1, jdata%block_rows
       CALL CopySparseMatrix(ABlocks(1,counter), TempABlocks(1,counter))
       CALL CopySparseMatrix(ABlocks(2,counter), TempABlocks(2,counter))
    END DO
    DO counter = 1, jdata%block_rows
       index = jdata%music_swap(counter)
       CALL CopySparseMatrix(TempABlocks(1,index), ABlocks(1,counter))
       CALL CopySparseMatrix(TempABlocks(2,index), ABlocks(2,counter))
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
          CALL SendMatrixSizes(SendLeft, jdata%send_left_partner, &
               & jdata%communicator, send_left_helper, &
               & jdata%send_left_tag)
          send_left_stage = send_left_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_left_helper)) THEN
             CALL SendMatrixData(SendLeft, jdata%send_left_partner, &
                  & jdata%communicator, send_left_helper, &
                  & jdata%send_left_tag)
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
          CALL RecvMatrixSizes(RecvLeft, jdata%recv_left_partner, &
               & jdata%communicator, recv_left_helper, &
               & jdata%recv_left_tag)
          recv_left_stage = recv_left_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_left_helper)) THEN
             CALL RecvMatrixData(RecvLeft, jdata%recv_left_partner, &
                  & jdata%communicator, recv_left_helper, &
                  & jdata%recv_left_tag)
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
          CALL SendMatrixSizes(SendRight, jdata%send_right_partner, &
               & jdata%communicator, send_right_helper, &
               & jdata%send_right_tag)
          send_right_stage = send_right_stage + 1
       CASE(1) !! Test Send Sizes
          IF (TestSendRecvSizeRequest(send_right_helper)) THEN
             CALL SendMatrixData(SendRight, jdata%send_right_partner, &
                  & jdata%communicator, send_right_helper, &
                  & jdata%send_right_tag)
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
          CALL RecvMatrixSizes(RecvRight, jdata%recv_right_partner, &
               & jdata%communicator, recv_right_helper, &
               & jdata%recv_right_tag)
          recv_right_stage = recv_right_stage + 1
       CASE(1) !! Test Receive Sizes
          IF (TestSendRecvSizeRequest(recv_right_helper)) THEN
             CALL RecvMatrixData(RecvRight, jdata%recv_right_partner, &
                  & jdata%communicator, recv_right_helper, &
                  & jdata%recv_right_tag)
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

    CALL SplitSparseMatrix(RecvLeft,jdata%block_rows,1,ABlocks(1,:))
    CALL SplitSparseMatrix(RecvRight,jdata%block_rows,1,ABlocks(2,:))

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

    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       CALL ComposeSparseMatrix(ABlocks(1:2,ind:ind+1), 2, 2, AMat)
       CALL Gemm(AMat, TargetV, TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(:,ind:ind+1))
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
    TYPE(SparseMatrix_t) :: TempMat
    TYPE(SparseMatrix_t) :: RecvMat
    TYPE(SparseMatrix_t) :: AMat
    INTEGER :: ind

    DO counter = 1, jdata%num_processes
       ind = (counter-1)*2 + 1
       IF (jdata%rank .EQ. counter-1) THEN
          CALL CopySparseMatrix(TargetV, RecvMat)
       END IF
       CALL BroadcastMatrix(RecvMat, jdata%communicator, counter - 1)

       CALL ComposeSparseMatrix(ABlocks(:,ind:ind+1), 2, 2, AMat)
       CALL Gemm(RecvMat, AMat, TempMat, threshold_in=threshold)
       CALL SplitSparseMatrix(TempMat, 2, 2, ABlocks(:,ind:ind+1))
       CALL DestructSparseMatrix(AMat)
    END DO

    CALL DestructSparseMatrix(RecvMat)
    CALL DestructSparseMatrix(TempMat)

  END SUBROUTINE ApplyToRows
END MODULE EigenSolversModule
