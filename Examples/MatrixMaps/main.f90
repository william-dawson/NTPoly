!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example which shows how to use the matrix mapping feature of NTPoly.
PROGRAM MatrixMapsProgram
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, ActivateLogger, &
       & DeactivateLogger, WriteElement, WriteHeader
  USE MatrixMapsModule, ONLY : MapMatrix_psr
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid, &
       & IsRoot
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructMatrixFromMatrixMarket, &
       & DestructMatrix, WriteMatrixToMatrixMarket
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  CHARACTER(len=80) :: input_matrix, output_matrix
  INTEGER :: process_slices
  !! Matrices
  TYPE(Matrix_ps) :: InMatrix, OutMatrix
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: II
  INTEGER :: provided, ierr
  INTEGER :: rank

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

  !! Process the input parameters.
  DO II = 1, COMMAND_ARGUMENT_COUNT(), 2
     CALL GET_COMMAND_ARGUMENT(II, argument)
     CALL GET_COMMAND_ARGUMENT(II+1, argument_value)
     SELECT CASE(argument)
     CASE('--input_matrix')
        input_matrix = argument_value
     CASE('--output_matrix')
        output_matrix = argument_value
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_slices)

  !! Print out parameters.
  IF (IsRoot()) THEN
     CALL ActivateLogger
  END IF
  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="input_matrix", VALUE=input_matrix)
  CALL WriteElement(key="output_matrix", VALUE=output_matrix)
  CALL WriteElement(key="process_slices", VALUE=process_slices)
  CALL ExitSubLog

  !! Read in the matrices from file.
  CALL ConstructMatrixFromMatrixMarket(InMatrix, input_matrix)

  !! Map
  CALL MapMatrix_psr(InMatrix, OutMatrix, TestFunction)

  !! Print the density matrix to file.
  CALL WriteMatrixToMatrixMarket(OutMatrix, output_matrix)

  !! Cleanup
  CALL DestructMatrix(InMatrix)
  CALL DestructMatrix(OutMatrix)

  !! Cleanup
  IF (IsRoot()) THEN
     CALL DeactivateLogger
  END IF
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This is the function we will map on to the matrix.
  FUNCTION TestFunction(row, column, val) RESULT(valid)
    INTEGER, INTENT(INOUT), OPTIONAL :: row
    INTEGER, INTENT(INOUT), OPTIONAL :: column
    REAL(NTREAL), INTENT(INOUT), OPTIONAL :: val
    LOGICAL :: valid

    IF (row .GE. column) THEN
       valid = .TRUE.
       val = val*2
    ELSE
       valid = .FALSE.
    END IF
  END FUNCTION TestFunction
END PROGRAM MatrixMapsProgram
