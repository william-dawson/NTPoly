!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example based on solving matrices based on premade files.
PROGRAM PremadeMatrixProgram
  USE DataTypesModule, ONLY : ntreal
  USE DensityMatrixSolversModule, ONLY : TRS2
  USE DistributedSparseMatrixModule, ONLY : ConstructFromMatrixMarket, &
       & WriteToMatrixMarket, DistributedSparseMatrix_t
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t
  USE LoggingModule
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation
  USE ProcessGridModule, ONLY : ConstructProcessGrid, IsRoot
  USE SquareRootSolversModule, ONLY : InverseSquareRoot
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  CHARACTER(len=80) :: hamiltonian_file
  CHARACTER(len=80) :: overlap_file
  CHARACTER(len=80) :: density_file_out
  INTEGER :: process_rows, process_columns, process_slices
  REAL(ntreal) :: threshold, convergence_threshold
  INTEGER :: number_of_electrons
  TYPE(IterativeSolverParameters_t) :: solver_parameters
  TYPE(Permutation_t) :: permutation
  !! Matrices
  TYPE(DistributedSparseMatrix_t) :: Hamiltonian, Overlap
  TYPE(DistributedSparseMatrix_t) :: ISQOverlap
  TYPE(DistributedSparseMatrix_t) :: Density
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: counter
  INTEGER :: provided, ierr
  INTEGER :: rank
  REAL(ntreal) :: chemical_potential

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

  !! Process the input parameters.
  DO counter=1,command_argument_count(),2
     CALL get_command_argument(counter,argument)
     CALL get_command_argument(counter+1,argument_value)
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
     CASE('--convergence_threshold')
        READ(argument_value,*) convergence_threshold
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="hamiltonian", text_value_in=hamiltonian_file)
  CALL WriteElement(key="overlap", text_value_in=overlap_file)
  CALL WriteElement(key="density", text_value_in=density_file_out)
  CALL WriteElement(key="process_rows", int_value_in=process_rows)
  CALL WriteElement(key="process_columns", int_value_in=process_columns)
  CALL WriteElement(key="process_slices", int_value_in=process_slices)
  CALL WriteElement(key="number_of_electrons", int_value_in=number_of_electrons)
  CALL WriteElement(key="threshold", float_value_in=threshold)
  CALL WriteElement(key="convergence_threshold", &
       & float_value_in=convergence_threshold)
  CALL ExitSubLog

  !! Read in the matrices from file.
  CALL ConstructFromMatrixMarket(Hamiltonian,hamiltonian_file)
  CALL ConstructFromMatrixMarket(Overlap,overlap_file)

  !! Set Up The Solver Parameters.
  CALL ConstructRandomPermutation(permutation, &
       & Hamiltonian%logical_matrix_dimension)
  solver_parameters = IterativeSolverParameters_t(&
       & converge_diff_in=convergence_threshold, threshold_in=threshold, &
       & BalancePermutation_in=permutation, be_verbose_in=.TRUE.)

  !! Call the solver routines.
  CALL InverseSquareRoot(Overlap, ISQOverlap, solver_parameters)
  CALL TRS2(Hamiltonian, ISQOverlap, number_of_electrons, &
       & Density, solver_parameters_in=solver_parameters, &
       & chemical_potential_out=chemical_potential)

  !! Print the density matrix to file.
  CALL WriteToMatrixMarket(Density,density_file_out)

  CALL MPI_Finalize(ierr)
END PROGRAM PremadeMatrixProgram
