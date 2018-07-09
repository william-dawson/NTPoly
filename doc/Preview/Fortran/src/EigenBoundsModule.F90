!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of a matrix's spectrum.
MODULE EigenBoundsModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, DotMatrix, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, GetMatrixTripletList, FillMatrixFromTripletList
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !> Uses Gershgorin's theorem.
  SUBROUTINE GershgorinBounds(this,min_value,max_value)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A lower bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: min_value
    !> An uppder bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !! Local Data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    !! Counters/Temporary
    INTEGER :: counter
    INTEGER :: local_column
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define triplet_list triplet_list_c
#include "SolverSupport/includes/GershgorinBounds.F90"
#undef triplet_list
    ELSE
#define triplet_list triplet_list_r
#include "SolverSupport/includes/GershgorinBounds.F90"
#undef triplet_list
    END IF
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !> Uses The Power Method.
  SUBROUTINE PowerBounds(this,max_value,solver_parameters_in)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> An upper bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Data
    TYPE(Matrix_ps) :: vector, vector2, TempMat
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: outer_counter
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
       solver_parameters%max_iterations = 10
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrix(vector, this)
    CALL ConstructEmptyMatrix(vector2, this)

    !! Guess Vector
    temp_list = TripletList_r()
    IF (this%process_grid%global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0
       CALL AppendToTripletList(temp_list,temp_triplet)
    END IF
    CALL FillMatrixFromTripletList(vector,temp_list)

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
       CALL MatrixMultiply(this,vector,vector2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !! x = x/||x||
       scale_value = 1.0/MatrixNorm(vector2)
       CALL ScaleMatrix(vector2,scale_value)

       !! Check if Converged
       CALL IncrementMatrix(vector2,vector,REAL(-1.0,NTREAL))
       norm_value = MatrixNorm(vector)

       CALL CopyMatrix(vector2,vector)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
    END IF

    !! Compute The Largest Eigenvalue
    CALL DotMatrix(vector, vector, scale_value)
    CALL MatrixMultiply(this,vector,vector2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL DotMatrix(vector, vector2, max_value)
    max_value = max_value / scale_value

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Max_Eigen_Value",float_value_in=max_value)
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(vector)
    CALL DestructMatrix(vector2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule
