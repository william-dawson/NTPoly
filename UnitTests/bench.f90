!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A minimal benchmark program
PROGRAM SimpleBench
  USE DataTypesModule
  USE LoadBalancerModule
  USE LoggingModule
  USE MatrixMapsModule
  USE PermutationModule
  USE ProcessGridModule
  USE PSMatrixAlgebraModule
  USE PSMatrixModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  INTEGER :: N
  REAL(NTREAL) :: thresh
  TYPE(Matrix_ps) :: A, B, D
  TYPE(Permutation_t) :: perm
  !! Temporary Variables
  INTEGER :: ierr, provided
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: II

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)

  !! Process the input parameters.
  DO II = 1, COMMAND_ARGUMENT_COUNT(), 2
     CALL GET_COMMAND_ARGUMENT(II, argument)
     CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
     SELECT CASE(argument)
     CASE('--N')
        READ(argument_value,*) N
     CASE('--thresh')
        READ(argument_value,*) thresh
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD)

  !! Print Out The Parameters
  IF (IsRoot()) THEN
     CALL ActivateLogger
  END IF
  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="N", VALUE=N)
  CALL WriteElement(key="thresh", VALUE=thresh)
  CALL ExitSubLog

  !! Timers
  CALL RegisterTimer("Construct")
  CALL RegisterTimer("Permute")
  CALL RegisterTimer("Multiply")

  !! Make the Matrices
  CALL StartTimer("Construct")
  CALL ConstructEmptyMatrix(A, N)
  CALL ConstructEmptyMatrix(B, N)
  CALL ConstructEmptyMatrix(D, N)
  CALL FillMatrixDense(D)
  CALL MapMatrix_psr(D, A, Shrink)
  CALL StopTimer("Construct")

  !! Permute
  CALL StartTimer("Permute")
  CALL ConstructRandomPermutation(perm, B%logical_matrix_dimension)
  CALL PermuteMatrix(A, A, perm)
  CALL StopTimer("Permute")

  !! Multiply
  CALL StartTimer("Multiply")
  CALL MatrixMultiply(A, A, B, threshold_in=thresh)
  CALL StopTimer("Multiply")
  
  !! Extra Data
  CALL WriteHeader("Matrix A")
  CALL EnterSubLog
  CALL PrintMatrixInformation(A)
  CALL ExitSubLog
  CALL WriteHeader("Matrix B")
  CALL EnterSubLog
  CALL PrintMatrixInformation(B)
  CALL ExitSubLog

  !! Cleanup
  CALL PrintAllTimers()
  CALL DestructMatrix(A)
  CALL DestructMatrix(B)
  CALL DestructMatrix(D)

  !! Cleanup
  IF (IsRoot()) THEN
     CALL DeactivateLogger
  END IF
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Shrink(row, column, val) RESULT(valid)
    INTEGER, INTENT(INOUT), OPTIONAL :: row
    INTEGER, INTENT(INOUT), OPTIONAL :: column
    REAL(NTREAL), INTENT(INOUT), OPTIONAL :: val
    LOGICAL :: valid
    REAL(NTREAL) :: r

    r = ABS(row - column + 0.5)
    val = -0.005 * EXP(-0.001 * r) / r
    valid = ABS(val) > thresh  * 10
  END FUNCTION Shrink
END PROGRAM SimpleBench

