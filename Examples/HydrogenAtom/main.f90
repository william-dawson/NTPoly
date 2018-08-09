!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example based on solving matrices based on a 1D hydrogen molecule.
PROGRAM HydrogenAtom
  USE DataTypesModule, ONLY : NTREAL
  USE DensityMatrixSolversModule, ONLY : TRS2
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  USE PSMatrixModule, ONLY : Matrix_ps, WriteMatrixToMatrixMarket, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, CopyMatrix, &
       & FillMatrixIdentity
  USE PSMatrixAlgebraModule, ONLY : IncrementMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t
  USE SquareRootSolversModule, ONLY : InverseSquareRoot
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  INTEGER :: grid_points
  CHARACTER(len=80) :: density_file_out
  INTEGER :: process_rows, process_columns, process_slices
  REAL(NTREAL) :: threshold, convergence_threshold
  TYPE(SolverParameters_t) :: solver_parameters
  !! MPI Variables
  INTEGER :: rank
  INTEGER :: total_processors
  INTEGER :: ierr
  INTEGER :: provided
  !! Linear Space
  INTEGER :: local_grid_points
  INTEGER, DIMENSION(:), ALLOCATABLE :: local_rows
  INTEGER :: start_row
  REAL(NTREAL), PARAMETER :: x_start = -6.28
  REAL(NTREAL), PARAMETER :: x_end = 6.28
  REAL(NTREAL) :: grid_spacing
  REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: x_values
  !! Matrices
  TYPE(Matrix_ps) :: KineticEnergy
  TYPE(Matrix_ps) :: PotentialEnergy
  TYPE(Matrix_ps) :: Hamiltonian
  TYPE(Matrix_ps) :: Identity
  TYPE(Matrix_ps) :: Density
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: counter

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, total_processors, ierr)

  !! Process the input parameters.
  DO counter=1,command_argument_count(),2
     CALL get_command_argument(counter,argument)
     CALL get_command_argument(counter+1,argument_value)
     SELECT CASE(argument)
     CASE('--convergence_threshold')
        READ(argument_value,*) convergence_threshold
     CASE('--grid_points')
        READ(argument_value,*) grid_points
     CASE('--density')
        density_file_out = argument_value
     CASE('--process_columns')
        READ(argument_value,*) process_columns
     CASE('--process_rows')
        READ(argument_value,*) process_rows
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     CASE('--threshold')
        READ(argument_value,*) threshold
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  !! Set Up The Solver Parameters.
  solver_parameters = SolverParameters_t( be_verbose_in=.TRUE., &
       & converge_diff_in=convergence_threshold, threshold_in=threshold)

  !! Divide The Work Amongst Processors.
  CALL DivideUpWork()

  !! Construct A Linear Space.
  CALL ConstructLinearSpace()

  !! Construct The Kinetic Energy Operator.
  CALL ConstructEmptyMatrix(KineticEnergy, grid_points)
  CALL FillKineticEnergy()

  !! Construct The Potential Energy Operator.
  CALL ConstructEmptyMatrix(PotentialEnergy, grid_points)
  CALL FillPotentialEnergy()

  !! Construct The Full Hamiltonian.
  CALL CopyMatrix(KineticEnergy, Hamiltonian)
  CALL IncrementMatrix(PotentialEnergy, Hamiltonian)

  !! Overlap Matrix is just the identity.
  CALL ConstructEmptyMatrix(Identity, Hamiltonian%actual_matrix_dimension)
  CALL FillMatrixIdentity(Identity)

  !! Call the solver routine.
  CALL TRS2(Hamiltonian, Identity, 2, Density, &
       & solver_parameters_in=solver_parameters)

  !! Print the density matrix to file.
  CALL WriteMatrixToMatrixMarket(Density,density_file_out)

  !! Cleanup
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DivideUpWork()
    local_grid_points = grid_points/total_processors
    start_row = local_grid_points * rank + 1
    !! Handles the edge case
    IF (rank .EQ. total_processors - 1) THEN
       local_grid_points = grid_points - rank*local_grid_points
    END IF
    ALLOCATE(local_rows(local_grid_points))
    fillrow: DO counter=1,local_grid_points
       local_rows(counter) = start_row + (counter-1)
    END DO fillrow
  END SUBROUTINE DivideUpWork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ConstructLinearSpace()
    !! Local Variables
    REAL(NTREAL) :: local_x_start

    ALLOCATE(x_values(local_grid_points))
    grid_spacing = (x_end - x_start)/(grid_points - 1)

    !! Fill in the x_values
    local_x_start = x_start + (start_row-1) * grid_spacing
    fill: DO counter=1,local_grid_points
       x_values(counter) = local_x_start + (counter-1)*grid_spacing
    END DO fill
  END SUBROUTINE ConstructLinearSpace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FillKineticEnergy()
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    INTEGER :: counter
    TYPE(Triplet_r) :: temp_value

    CALL ConstructTripletList(triplet_list)

    !! Fill The Triplet List.
    fill: DO counter = 1, local_grid_points
       !! Stencil point 1.
       temp_value%index_row = start_row + counter - 1
       IF (temp_value%index_row .GT. 2) THEN
          temp_value%index_column = temp_value%index_row - 2
          temp_value%point_value  = (-0.5)*(-1.0/(12.0*grid_spacing**2))
          CALL AppendToTripletList(triplet_list, temp_value)
       END IF
       !! Stencil point 2.
       IF (temp_value%index_row .GT. 1) THEN
          temp_value%index_column = temp_value%index_row - 1
          temp_value%point_value  = (-0.5)*(16.0/(12.0*grid_spacing**2))
          CALL AppendToTripletList(triplet_list, temp_value)
       END IF
       !! Stencil point 3.
       temp_value%index_column = temp_value%index_row
       temp_value%point_value  = (-0.5)*(-30.0/(12.0*grid_spacing**2))
       CALL AppendToTripletList(triplet_list, temp_value)
       !! Stencil point 4.
       IF (temp_value%index_row .LT. grid_points) THEN
          temp_value%index_column = temp_value%index_row + 1
          temp_value%point_value  = (-0.5)*(16.0/(12.0*grid_spacing**2))
          CALL AppendToTripletList(triplet_list, temp_value)
       END IF
       !! Stencil point 5.
       IF (temp_value%index_row .LT. grid_points - 1) THEN
          temp_value%index_column = temp_value%index_row + 2
          temp_value%point_value  = (-0.5)*(-1.0/(12.0*grid_spacing**2))
          CALL AppendToTripletList(triplet_list, temp_value)
       END IF
    END DO fill
    CALL FillMatrixFromTripletList(KineticEnergy, triplet_list)
  END SUBROUTINE FillKineticEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FillPotentialEnergy()
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    INTEGER :: counter
    TYPE(Triplet_r) :: temp_value

    CALL ConstructTripletList(triplet_list)
    fill: DO counter = 1, local_grid_points
       temp_value%index_row = start_row + counter - 1
       temp_value%index_column = start_row + counter - 1
       temp_value%point_value = -1.0/ABS(x_values(counter))
       CALL AppendToTripletList(triplet_list, temp_value)
    END DO fill
    CALL FillMatrixFromTripletList(PotentialEnergy, triplet_list)
  END SUBROUTINE FillPotentialEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM HydrogenAtom
