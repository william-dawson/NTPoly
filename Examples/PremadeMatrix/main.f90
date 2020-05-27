!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example based on solving matrices based on premade files.
PROGRAM PremadeMatrixProgram
  USE DataTypesModule, ONLY : NTREAL
  USE DensityMatrixSolversModule, ONLY : TRS2
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, WriteHeader
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation, &
       & DestructPermutation
  USE ProcessGridModule, ONLY : ConstructProcessGrid, IsRoot, &
       & DestructProcessGrid
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructMatrixFromMatrixMarket, &
       & WriteMatrixToMatrixMarket, DestructMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t
  USE SquareRootSolversModule, ONLY : InverseSquareRoot
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  CHARACTER(len=80) :: hamiltonian_file
  CHARACTER(len=80) :: overlap_file
  CHARACTER(len=80) :: density_file_out
  INTEGER :: process_rows, process_columns, process_slices
  REAL(NTREAL) :: threshold
  REAL(NTREAL) :: converge_overlap, converge_density
  INTEGER :: number_of_electrons
  TYPE(SolverParameters_t) :: solver_parameters
  TYPE(Permutation_t) :: permutation
  !! Matrices
  TYPE(Matrix_ps) :: Hamiltonian, Overlap
  TYPE(Matrix_ps) :: ISQOverlap
  TYPE(Matrix_ps) :: Density
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: counter
  INTEGER :: provided, ierr
  INTEGER :: rank
  REAL(NTREAL) :: chemical_potential

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

  !! Process the input parameters.
  DO counter=1,COMMAND_ARGUMENT_COUNT(),2
     CALL GET_COMMAND_ARGUMENT(counter,argument)
     CALL GET_COMMAND_ARGUMENT(counter+1,argument_value)
     SELECT CASE(argument)
     CASE('--hamiltonian')
        hamiltonian_file = argument_value
     CASE('--overlap')
        overlap_file = argument_value
     CASE('--density')
        density_file_out = argument_value
     CASE('--process_rows')
        READ(argument_value,*) process_rows
     CASE('--process_columns')
        READ(argument_value,*) process_columns
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     CASE('--number_of_electrons')
        READ(argument_value,*) number_of_electrons
     CASE('--threshold')
        READ(argument_value,*) threshold
     CASE('--converge_density')
        READ(argument_value,*) converge_density
     CASE('--converge_overlap')
        READ(argument_value,*) converge_overlap
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="hamiltonian", VALUE=hamiltonian_file)
  CALL WriteElement(key="overlap", VALUE=overlap_file)
  CALL WriteElement(key="density", VALUE=density_file_out)
  CALL WriteElement(key="process_rows", VALUE=process_rows)
  CALL WriteElement(key="process_columns", VALUE=process_columns)
  CALL WriteElement(key="process_slices", VALUE=process_slices)
  CALL WriteElement(key="number_of_electrons", VALUE=number_of_electrons)
  CALL WriteElement(key="threshold", VALUE=threshold)
  CALL WriteElement(key="converge_overlap", VALUE=converge_overlap)
  CALL WriteElement(key="converge_density", VALUE=converge_density)
  CALL ExitSubLog

  !! Read in the matrices from file.
  CALL ConstructMatrixFromMatrixMarket(Hamiltonian,hamiltonian_file)
  CALL ConstructMatrixFromMatrixMarket(Overlap,overlap_file)

  !! Set Up The Solver Parameters.
  CALL ConstructRandomPermutation(permutation, &
       & Hamiltonian%logical_matrix_dimension)
  solver_parameters = SolverParameters_t(&
       & converge_diff_in=converge_overlap, threshold_in=threshold, &
       & BalancePermutation_in=permutation, be_verbose_in=.TRUE.)

  !! Call the solver routines.
  CALL InverseSquareRoot(Overlap, ISQOverlap, solver_parameters)

  !! Change the solver variable for computing the density matrix.
  solver_parameters%converge_diff=converge_density

  !! Compute the density matrix.
  CALL TRS2(Hamiltonian, ISQOverlap, number_of_electrons, &
       & Density, solver_parameters_in=solver_parameters, &
       & chemical_potential_out=chemical_potential)

  !! Print the density matrix to file.
  CALL WriteMatrixToMatrixMarket(Density,density_file_out)

  !! Cleanup
  CALL DestructPermutation(permutation)
  CALL DestructMatrix(Overlap)
  CALL DestructMatrix(ISQOverlap)
  CALL DestructMatrix(Hamiltonian)
  CALL DestructMatrix(Density)

  !! Cleanup
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
END PROGRAM PremadeMatrixProgram
