!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Study the exponential of a nonsymmetric graph using the technique of
!! Guo, Krystal, and Bojan Mohar. "Hermitian adjacency matrix of digraphs and
!! mixed graphs." Journal of Graph Theory 85, no. 1 (2017): 217-248.
PROGRAM ComplexMatrix
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE ExponentialSolversModule, ONLY : ComputeExponential
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, WriteHeader
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructMatrixFromMatrixMarket, &
       & DestructMatrix, TransposeMatrix, GetMatrixTripletList, &
       & FillMatrixFromTripletList, ConstructEmptyMatrix, ConjugateMatrix, &
       & PrintMatrix, CopyMatrix, WriteMatrixToMatrixMarket
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  USE SolverParametersModule, ONLY : SolverParameters_t
  USE TripletListModule, ONLY : TripletList_r, DestructTripletList, &
       & GetTripletAt, TripletList_c, ConstructTripletList, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  INTEGER :: process_rows, process_columns, process_slices
  CHARACTER(len=80) :: input_file
  CHARACTER(len=80) :: exponential_file
  REAL(NTREAL) :: threshold
  TYPE(SolverParameters_t) :: solver_parameters
  !! Matrices
  TYPE(Matrix_ps) :: InMat
  TYPE(Matrix_ps) :: GMat
  TYPE(Matrix_ps) :: ExMat
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: counter
  INTEGER :: provided, ierr

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)

  !! Process The Input
  !! Process the input parameters.
  DO counter=1,COMMAND_ARGUMENT_COUNT(),2
     CALL GET_COMMAND_ARGUMENT(counter,argument)
     CALL GET_COMMAND_ARGUMENT(counter+1,argument_value)
     SELECT CASE(argument)
     CASE('--input_file')
        input_file = argument_value
     CASE('--exponential_file')
        exponential_file = argument_value
     CASE('--process_rows')
        READ(argument_value,*) process_rows
     CASE('--process_columns')
        READ(argument_value,*) process_columns
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     CASE('--threshold')
        READ(argument_value,*) threshold
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="input_file", VALUE=input_file)
  CALL WriteElement(key="exponential_file", VALUE=exponential_file)
  CALL WriteElement(key="process_rows", VALUE=process_rows)
  CALL WriteElement(key="process_columns", VALUE=process_columns)
  CALL WriteElement(key="process_slices", VALUE=process_slices)
  CALL WriteElement(key="threshold", VALUE=threshold)
  CALL ExitSubLog

  !! Construct The Hermitian Matrix
  CALL ConstructMatrixFromMatrixMarket(InMat, input_file)
  CALL ConstructGuoMatrix(InMat, GMat)
  CALL ScaleMatrix(GMat, 0.5_NTREAL)

  !! Compute The Exponential
  solver_parameters = SolverParameters_t(threshold_in=threshold)
  CALL ComputeExponential(GMat, ExMat, solver_parameters)

  !! Write To File
  CALL WriteMatrixToMatrixMarket(ExMat, exponential_file)

  !! Cleanup
  CALL DestructMatrix(InMat)
  CALL DestructMatrix(GMat)
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)

CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the Hermitian matrix from the nonsymmetric matrix
  SUBROUTINE ConstructGuoMatrix(Input, Output)
    !> Input matrix
    TYPE(Matrix_ps), INTENT(IN) :: Input
    !> The hermitian output matrix
    TYPE(Matrix_ps), INTENT(INOUT) :: Output
    !! Local Variables
    TYPE(Matrix_ps) :: Guide
    TYPE(Matrix_ps) :: SMatrix
    TYPE(Matrix_ps) :: CMatrix
    TYPE(TripletList_r) :: tlist
    TYPE(TripletList_c) :: clist
    !! Temporary Variables
    TYPE(Triplet_r) :: temp
    TYPE(Triplet_c) :: temp_c
    INTEGER :: II

    !! Symmetric The Input Matrix
    CALL GetMatrixTripletList(Input, tlist)
    CALL SymmetrizeTripletList(tlist, MM_SYMMETRIC)
    CALL ConstructEmptyMatrix(SMatrix, Input)
    CALL FillMatrixFromTripletList(SMatrix, tlist)

    !! Construct The Guide Matrix
    CALL CopyMatrix(SMatrix, Guide)
    CALL IncrementMatrix(Input, Guide, alpha_in=-1.0_NTREAL)

    !! Now iterate over the entries in the guide matrix.
    CALL GetMatrixTripletList(Guide, tlist)
    CALL ConstructTripletList(clist)
    DO II = 1, tlist%CurrentSize
       CALL GetTripletAt(tlist, II, temp)
       temp_c%index_row = temp%index_row
       temp_c%index_column = temp%index_column
       temp_c%point_value = (0.0,1.0)
       CALL AppendToTripletList(clist, temp_c)
    END DO
    CALL ConstructEmptyMatrix(CMatrix, Input)
    CALL FillMatrixFromTripletList(CMatrix, clist, preduplicated_in=.TRUE.)

    !! Now add it all together
    CALL TransposeMatrix(CMatrix, Output)
    CALL ConjugateMatrix(Output)
    CALL IncrementMatrix(CMatrix, Output)
    CALL IncrementMatrix(SMatrix, Output)

    !! Cleanup
    CALL DestructMatrix(SMatrix)
    CALL DestructMatrix(CMatrix)
    CALL DestructMatrix(Guide)
    CALL DestructTripletList(clist)
    CALL DestructTripletList(tlist)
  END SUBROUTINE ConstructGuoMatrix
END PROGRAM ComplexMatrix
