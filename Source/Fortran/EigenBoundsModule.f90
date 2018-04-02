!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of a matrix's spectrum.
MODULE EigenBoundsModule
  USE DataTypesModule
  USE DenseMatrixModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE FixedSolversModule
  USE IterativeSolversModule
  USE LoggingModule
  USE MatrixGatherModule
  USE ProcessGridModule
  USE SparseMatrixModule
  USE SparseMatrixAlgebraModule
  USE TripletListModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
  PUBLIC :: DistributedEigenDecomposition
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !! Uses Gershgorin's theorem.
  !! @param[in] this the matrix to compute the min/max of.
  !! @param[out] min_value a lower bound on the eigenspectrum.
  !! @param[out] max_value an uppder bound on the eigenspectrum.
  SUBROUTINE GershgorinBounds(this,min_value,max_value)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    REAL(NTREAL), INTENT(OUT) :: min_value, max_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    TYPE(TripletList_t) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    INTEGER :: local_column

    !! Allocate Space For Result
    ALLOCATE(per_column_min(this%local_columns))
    ALLOCATE(per_column_max(this%local_columns))

    !! Compute The Local Contribution
    per_column_min = 0
    per_column_max = 0
    CALL GetTripletList(this, triplet_list)
    DO counter = 1, triplet_list%CurrentSize
       local_column = triplet_list%data(counter)%index_column - &
            & this%start_column + 1
       IF (triplet_list%data(counter)%index_row .EQ. &
            & triplet_list%data(counter)%index_column) THEN
          per_column_min(local_column) = per_column_min(local_column) + &
               & triplet_list%data(counter)%point_value
          per_column_max(local_column) = per_column_max(local_column) + &
               & triplet_list%data(counter)%point_value
       ELSE
          per_column_min(local_column) = per_column_min(local_column) - &
               & ABS(triplet_list%data(counter)%point_value)
          per_column_max(local_column) = per_column_max(local_column) + &
               & ABS(triplet_list%data(counter)%point_value)
       END IF
    END DO

    !! Sum Along Columns
    CALL MPI_Allreduce(MPI_IN_PLACE,per_column_min,SIZE(per_column_min), &
         & MPINTREAL,MPI_SUM,column_comm,grid_error)
    CALL MPI_Allreduce(MPI_IN_PLACE,per_column_max,SIZE(per_column_max), &
         & MPINTREAL,MPI_SUM,column_comm,grid_error)

    min_value = MINVAL(per_column_min)
    max_value = MAXVAL(per_column_max)

    CALL MPI_Allreduce(MPI_IN_PLACE,min_value,1,MPINTREAL,MPI_MIN, &
         & row_comm, grid_error)
    CALL MPI_Allreduce(MPI_IN_PLACE,max_value,1,MPINTREAL,MPI_MAX, &
         & row_comm, grid_error)

    CALL DestructTripletList(triplet_list)
    DEALLOCATE(per_column_min)
    DEALLOCATE(per_column_max)
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !! Uses The Power Method.
  !! @param[in] this the matrix to compute the min/max of.
  !! @param[out] max_value an upper bound on the eigenspectrum.
  !! @param[inout] solver_parameters_in solver parameters (optional).
  SUBROUTINE PowerBounds(this,max_value,solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    REAL(NTREAL), INTENT(OUT) :: max_value
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Data
    TYPE(DistributedSparseMatrix_t) :: vector, vector2, TempMat
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_t) :: temp_list
    TYPE(Triplet_t) :: temp_triplet
    INTEGER :: outer_counter
    TYPE(DistributedMatrixMemoryPool_t) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
       solver_parameters%max_iterations = 10
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyDistributedSparseMatrix(vector, &
         & this%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(vector2, &
         & this%actual_matrix_dimension)

    !! Guess Vector
    CALL ConstructTripletList(temp_list)
    IF (global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0
       CALL AppendToTripletList(temp_list,temp_triplet)
    END IF
    CALL FillFromTripletList(vector,temp_list)

    !! Iterate
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       !! x = Ax
       CALL DistributedGemm(this,vector,vector2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !! x = x/||x||
       scale_value = 1.0/DistributedSparseNorm(vector2)
       CALL ScaleDistributedSparseMatrix(vector2,scale_value)

       !! Check if Converged
       CALL IncrementDistributedSparseMatrix(vector2,vector,REAL(-1.0,NTREAL))
       norm_value = DistributedSparseNorm(vector)

       CALL CopyDistributedSparseMatrix(vector2,vector)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
    END IF

    !! Compute The Largest Eigenvalue
    scale_value = DotDistributedSparseMatrix(vector,vector)
    CALL DistributedGemm(this,vector,vector2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    max_value = DotDistributedSparseMatrix(vector,vector2)
    max_value = max_value / scale_value

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Max_Eigen_Value",float_value_in=max_value)
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructDistributedSparseMatrix(vector)
    CALL DestructDistributedSparseMatrix(vector2)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a matrix using Jacobi's method.
  !! @param[in] this the matrix to compute the eigenvectors of.
  !! @param[out] eigenvectors the eigenvectors of the matrix
  !! @param[inout] solver_parameters_in solver parameters (optional).
  SUBROUTINE DistributedEigenDecomposition(this, eigenvectors, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: this
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: eigenvectors
    TYPE(FixedSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(FixedSolverParameters_t) :: solver_parameters
    !! Local Variables
    TYPE(SparseMatrix_t) :: local_a, local_v, local_p
    TYPE(SparseMatrix_t) :: temp, temp2
    TYPE(SparseMatrix_t), DIMENSION(4) :: send_block, recv_block
    !! Temporary Variables
    INTEGER :: iteration
    INTEGER :: block_counter
    INTEGER :: root_row, root_column

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = FixedSolverParameters_t()
    END IF

    !! Setup
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, &
         & this%actual_matrix_dimension)
    CALL FillDistributedIdentity(eigenvectors)
    CALL MergeLocalBlocks(eigenvectors, local_v)
    CALL MergeLocalBlocks(this, local_a)
    root_row = my_row
    root_column = my_column

    !! Main Loop
    DO iteration = 1, 2
       DO block_counter = 1, 1
          !! Compute the eigenvectors of a local block
          IF (my_row .EQ. my_column) THEN
             CALL DenseEigenDecomposition(local_a, local_p, &
                  & solver_parameters%threshold)
          END IF

          !! Broadcast P
          CALL BroadcastMatrix(local_p, row_comm, root_row)
          CALL BroadcastMatrix(local_p, column_comm, root_column)

          !! Multiply Blocks
          CALL Gemm(local_a, local_p, temp, &
               & threshold_in=solver_parameters%threshold)
          CALL TransposeSparseMatrix(local_p, temp2)
          CALL Gemm(temp2, temp, local_a, &
               & threshold_in=solver_parameters%threshold)
          CALL Gemm(local_v, local_p, temp, &
               & threshold_in=solver_parameters%threshold)
          CALL CopySparseMatrix(temp, local_v)

          !! Exchange Blocks
          ! CALL SplitFour(local_a, send_block)
          ! IF (my_row .EQ. my_column .AND. my_row .EQ. 1) THEN
          !   CALL CopyDistributedSparseMatrix(send_block(1),recv_block(1))
          ! END IF
          ! CALL MergeFour(recv_block, local_a)
       END DO
    END DO

    !! Compute Full Matrix
    CALL ConstructEmptyDistributedSparseMatrix(eigenvectors, this%actual_matrix_dimension)
    CALL SplitToLocalBlocks(eigenvectors, local_v)

    !! Cleanup
    CALL DestructSparseMatrix(temp)
    CALL DestructSparseMatrix(temp2)
    CALL DestructSparseMatrix(local_v)
    CALL DestructSparseMatrix(local_a)
    CALL DestructSparseMatrix(local_p)
    CALL DestructSparseMatrix(send_block(1))
    CALL DestructSparseMatrix(send_block(2))
    CALL DestructSparseMatrix(send_block(3))
    CALL DestructSparseMatrix(send_block(4))
    CALL DestructSparseMatrix(recv_block(1))
    CALL DestructSparseMatrix(recv_block(2))
    CALL DestructSparseMatrix(recv_block(3))
    CALL DestructSparseMatrix(recv_block(4))

  END SUBROUTINE DistributedEigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE SplitFour(matrix, matrix_array)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: matrix
    TYPE(SparseMatrix_t), DIMENSION(4), INTENT(INOUT) :: matrix_array
    !! Local
    TYPE(SparseMatrix_t) :: temp
    TYPE(SparseMatrix_t), DIMENSION(2) :: temp_rows
    TYPE(SparseMatrix_t), DIMENSION(2) :: temp_columns

    CALL SplitSparseMatrixColumns(matrix, 2, temp_columns)

    CALL TransposeSparseMatrix(temp_columns(1), temp)
    CALL SplitSparseMatrixColumns(temp, 2, temp_rows)
    CALL TransposeSparseMatrix(temp_rows(1), matrix_array(1))
    CALL TransposeSparseMatrix(temp_rows(2), matrix_array(3))

    CALL TransposeSparseMatrix(temp_columns(2), temp)
    CALL SplitSparseMatrixColumns(temp, 2, temp_rows)
    CALL TransposeSparseMatrix(temp_rows(1), matrix_array(2))
    CALL TransposeSparseMatrix(temp_rows(2), matrix_array(4))

    !! Cleanup
    CALL DestructSparseMatrix(temp_rows(1))
    CALL DestructSparseMatrix(temp_rows(2))
    CALL DestructSparseMatrix(temp_columns(1))
    CALL DestructSparseMatrix(temp_columns(2))
  END SUBROUTINE SplitFour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE MergeFour(matrix_array, matrix)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(4), INTENT(IN) :: matrix_array
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    !! Local
    TYPE(SparseMatrix_t) :: temp
    TYPE(SparseMatrix_t), DIMENSION(2) :: temp_rows

    CALL ComposeSparseMatrixColumns(matrix_array(1:2), temp)
    CALL TransposeSparseMatrix(temp, temp_rows(1))

    CALL ComposeSparseMatrixColumns(matrix_array(3:4), temp)
    CALL TransposeSparseMatrix(temp, temp_rows(2))

    CALL ComposeSparseMatrixColumns(temp_rows, temp)
    CALL TransposeSparseMatrix(temp, matrix)

    !! Cleanup
    CALL DestructSparseMatrix(temp)
    CALL DestructSparseMatrix(temp_rows(1))
    CALL DestructSparseMatrix(temp_rows(2))
  END SUBROUTINE MergeFour
END MODULE EigenBoundsModule
