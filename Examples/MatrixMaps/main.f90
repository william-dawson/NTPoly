!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example which shows how to use the matrix mapping feature of NTPoly.
PROGRAM MatrixMapsProgram
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, WriteHeader
  USE MatrixMapsModule, ONLY : MapMatrix_psr
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
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
  INTEGER :: counter
  INTEGER :: provided, ierr
  INTEGER :: rank

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

  !! Process the input parameters.
  DO counter=1,command_argument_count(),2
     CALL get_command_argument(counter,argument)
     CALL get_command_argument(counter+1,argument_value)
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
  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="input_matrix", value=input_matrix)
  CALL WriteElement(key="output_matrix", value=output_matrix)
  CALL WriteElement(key="process_slices", value=process_slices)
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
