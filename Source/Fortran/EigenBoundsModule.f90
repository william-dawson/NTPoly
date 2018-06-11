!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of a matrix's spectrum.
MODULE EigenBoundsModule
  USE DataTypesModule
  USE MatrixMemoryPoolDModule
  USE MatrixDSAlgebraModule
  USE MatrixDSModule
  USE IterativeSolversModule
  USE LoggingModule
  USE TripletListRModule
  USE TripletRModule
  USE TripletCModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !! Uses Gershgorin's theorem.
  !! @param[in] this the matrix to compute the min/max of.
  !! @param[out] min_value a lower bound on the eigenspectrum.
  !! @param[out] max_value an uppder bound on the eigenspectrum.
  SUBROUTINE GershgorinBounds(this,min_value,max_value)
    !! Parameters
    TYPE(Matrix_ds), INTENT(IN) :: this
    REAL(NTREAL), INTENT(OUT) :: min_value, max_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    TYPE(TripletList_r) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    INTEGER :: local_column
    INTEGER :: ierr

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
         & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE,per_column_max,SIZE(per_column_max), &
         & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)

    min_value = MINVAL(per_column_min)
    max_value = MAXVAL(per_column_max)

    CALL MPI_Allreduce(MPI_IN_PLACE,min_value,1,MPINTREAL,MPI_MIN, &
         & this%process_grid%row_comm, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE,max_value,1,MPINTREAL,MPI_MAX, &
         & this%process_grid%row_comm, ierr)

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
    TYPE(Matrix_ds), INTENT(IN) :: this
    REAL(NTREAL), INTENT(OUT) :: max_value
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Data
    TYPE(Matrix_ds) :: vector, vector2, TempMat
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: outer_counter
    TYPE(MatrixMemoryPool_d) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
       solver_parameters%max_iterations = 10
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrixDS(vector, &
         & this%actual_matrix_dimension, this%process_grid)
    CALL ConstructEmptyMatrixDS(vector2, &
         & this%actual_matrix_dimension, this%process_grid)

    !! Guess Vector
    CALL ConstructTripletList(temp_list)
    IF (this%process_grid%global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0
       CALL AppendToTripletList(temp_list,temp_triplet)
    END IF
    CALL FillFromTripletList(vector,temp_list)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
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

       CALL CopyMatrixDS(vector2,vector)

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
    CALL DestructMatrixDS(vector)
    CALL DestructMatrixDS(vector2)
    CALL DestructMatrixDS(TempMat)
    CALL DestructMatrixMemoryPoolD(pool)
  END SUBROUTINE PowerBounds
END MODULE EigenBoundsModule
