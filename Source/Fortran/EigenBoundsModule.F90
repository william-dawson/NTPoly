!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of the spectrum of a matrix.
MODULE EigenBoundsModule
  USE ConvergenceMonitor, ONLY : ConstructMonitor, CheckConverged, AppendValue
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, DotMatrix, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, GetMatrixTripletList, FillMatrixFromTripletList
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, DestructTripletList, ConstructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !> Uses the Gershgorin theorem.
  SUBROUTINE GershgorinBounds(this, min_value, max_value)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A lower bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: min_value
    !> An uppder bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !! Local Data
    TYPE(TripletList_r) :: tlist_r
    TYPE(TripletList_c) :: tlist_c
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    !! Counters/Temporary
    INTEGER :: II
    INTEGER :: local_column
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define tlist tlist_c
#include "solver_includes/GershgorinBounds.f90"
#undef tlist
    ELSE
#define tlist tlist_r
#include "solver_includes/GershgorinBounds.f90"
#undef tlist
    END IF
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !> Uses The Power Method.
  SUBROUTINE PowerBounds(this, max_value, solver_parameters_in)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> An upper bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Data
    REAL(NTREAL), DIMENSION(3) :: ritz_values
    REAL(NTREAL), DIMENSION(3) :: aitken_values
    REAL(NTREAL) :: num, den
    TYPE(Matrix_ps) :: vector, vector2
    REAL(NTREAL) :: scale_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: II
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
       params%max_iterations = 10
    END IF
    CALL ConstructMonitor(params%monitor, &
         & automatic_in = params%monitor_convergence, &
         & tight_cutoff_in=params%converge_diff)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintParameters(params)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrix(vector, this)
    CALL ConstructEmptyMatrix(vector2, this)

    !! Guess Vector
    CALL ConstructTripletList(temp_list)
    IF (this%start_row .EQ. 1) THEN
       DO II = this%start_column, this%end_column - 1
          temp_triplet%index_row = 1
          temp_triplet%index_column = II
          temp_triplet%point_value = &
               & 1.0_NTREAL / vector%actual_matrix_dimension
          CALL AppendToTripletList(temp_list, temp_triplet)
       END DO
    END IF
    CALL FillMatrixFromTripletList(vector, temp_list, &
         & preduplicated_in=.TRUE., prepartitioned_in=.TRUE.)

    !! Iterate
    ritz_values = 0
    aitken_values = 0
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    DO II = 1, params%max_iterations
       !! x = Ax
       CALL MatrixMultiply(this, vector, vector2, &
            & threshold_in = params%threshold, memory_pool_in = pool)

       !! Compute the ritz value
       CALL DotMatrix(vector, vector, scale_value)
       CALL DotMatrix(vector, vector2, max_value)
       max_value = max_value / scale_value

       !! x = x/||x||
       scale_value = 1.0 / MatrixNorm(vector2)
       CALL ScaleMatrix(vector2, scale_value)
       CALL CopyMatrix(vector2, vector)

       !! Aitken Extrapolation
       ritz_values(1) = ritz_values(2)
       ritz_values(2) = ritz_values(3)
       ritz_values(3) = max_value
       aitken_values(1) = aitken_values(2)
       aitken_values(2) = aitken_values(3)
       IF (II .GE. 3) THEN
          num = ritz_values(3)*ritz_values(1) - ritz_values(2)**2
          den = ritz_values(3) - 2 * ritz_values(2) + ritz_values(1)
          !! Avoid division by zero
          IF (ABS(den) .GT. 1E-14_NTREAL) THEN
             aitken_values(3) = num / den
          ELSE
             aitken_values(3) = ritz_values(3)
          END IF
       ELSE
          aitken_values(3) = ritz_values(3)
       END IF

       !! Check if Converged - pass the negative value because we are looking
       !! for the largest eigenvalue value.
       CALL AppendValue(params%monitor,  &
            & - (aitken_values(3) - aitken_values(2)))
       IF (CheckConverged(params%monitor, params%be_verbose)) THEN
          !! Make sure the two estimates vaguely agree
          IF(ABS(aitken_values(3) - ritz_values(3)) .LT. &
               & params%monitor%loose_cutoff) THEN
             EXIT
          END IF
       END IF
       IF (params%be_verbose) THEN
          CALL EnterSubLog
          CALL WriteElement(key="Estimate", VALUE=ritz_values(3))
          CALL WriteElement(key="Aitken Estimate", VALUE=aitken_values(3))
          CALL ExitSubLog
       END IF

    END DO
    max_value = aitken_values(3)
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II - 1)
       CALL WriteElement(key = "Max Eigen Value", VALUE = aitken_values(3))
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(vector)
    CALL DestructMatrix(vector2)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule
