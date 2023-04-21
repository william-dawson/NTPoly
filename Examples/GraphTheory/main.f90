!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> An example based on graph theory.
PROGRAM GraphTheory
  USE DataTypesModule, ONLY : NTREAL
  USE InverseSolversModule, ONLY : Invert
  USE LoggingModule, ONLY : ActivateLogger, DeactivateLogger, WriteElement, &
       & WriteHeader, EnterSubLog, ExitSubLog
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid, &
       & IsRoot
  USE PSMatrixModule, ONLY : Matrix_ps, WriteMatrixToMatrixMarket, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, DestructMatrix, &
       & CopyMatrix, FillMatrixIdentity
  USE PSMatrixAlgebraModule, ONLY: IncrementMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & ConstructSolverParameters, DestructSolverParameters
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & SetTripletAt, AppendToTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  INTEGER :: process_rows, process_columns, process_slices
  INTEGER :: number_of_nodes
  INTEGER :: extra_connections
  CHARACTER(len=80) :: output_file
  REAL(NTREAL) :: threshold
  REAL(NTREAL) :: convergence_threshold
  REAL(NTREAL) :: attenuation
  TYPE(SolverParameters_t) :: solver_parameters
  !! MPI Variables
  INTEGER :: rank
  INTEGER :: total_processors
  INTEGER :: ierr
  INTEGER :: provided
  !! Matrices
  TYPE(Matrix_ps) :: NetworkMat
  TYPE(Matrix_ps) :: ResultMat
  !! Local Part of the Matrix
  INTEGER :: number_of_local_nodes
  INTEGER, DIMENSION(:), ALLOCATABLE :: local_nodes
  INTEGER :: starting_node
  INTEGER :: ending_node
  !! Temporary Values
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: II

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, total_processors, ierr)

  !! Process the input parameters.
  DO II = 1, COMMAND_ARGUMENT_COUNT(), 2
     CALL GET_COMMAND_ARGUMENT(II, argument)
     CALL GET_COMMAND_ARGUMENT(II + 1, argument_value)
     SELECT CASE(argument)
     CASE('--output_file')
        output_file = argument_value
     CASE('--process_rows')
        READ(argument_value,*) process_rows
     CASE('--process_columns')
        READ(argument_value,*) process_columns
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     CASE('--number_of_nodes')
        READ(argument_value,*) number_of_nodes
     CASE('--extra_connections')
        READ(argument_value,*) extra_connections
     CASE('--threshold')
        READ(argument_value,*) threshold
     CASE('--convergence_threshold')
        READ(argument_value,*) convergence_threshold
     CASE('--attenuation')
        READ(argument_value,*) attenuation
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  !! Write Out Parameters
  IF (IsRoot()) THEN
     CALL ActivateLogger
  END IF
  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="number_of_nodes", VALUE=number_of_nodes)
  CALL WriteElement(key="extra_connections", VALUE=extra_connections)
  CALL WriteElement(key="threshold", VALUE=threshold)
  CALL WriteElement(key="convergence_threshold", VALUE=convergence_threshold)
  CALL WriteElement(key="attenuation", VALUE=attenuation)
  CALL WriteElement(key="process_rows", VALUE=process_rows)
  CALL WriteElement(key="process_columns", VALUE=process_columns)
  CALL WriteElement(key="process_slices", VALUE=process_slices)
  CALL WriteElement(key="output_file", VALUE=output_file)
  CALL ExitSubLog

  !! Set Up The Solver Parameters.
  CALL ConstructSolverParameters(solver_parameters, be_verbose_in=.TRUE., &
       & converge_diff_in=convergence_threshold, threshold_in=threshold)

  CALL DivideUpWork

  !! Fill The Matrix
  CALL ConstructEmptyMatrix(NetworkMat, number_of_nodes)
  CALL ConstructEmptyMatrix(ResultMat, number_of_nodes)
  CALL FillMatrix

  !! Solve
  CALL SolveMatrix

  !! Print the density matrix to file.
  CALL WriteMatrixToMatrixMarket(ResultMat,output_file)

  !! Cleanup
  CALL DestructMatrix(NetworkMat)
  CALL DestructMatrix(ResultMat)

  !! Cleanup
  IF (IsRoot()) THEN
     CALL DeactivateLogger
  END IF
  CALL DestructSolverParameters(solver_parameters)
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DivideUpWork()
    number_of_local_nodes = number_of_nodes/total_processors
    starting_node = number_of_local_nodes * rank + 1
    !! Handles the edge case
    IF (rank .EQ. total_processors - 1) THEN
       number_of_local_nodes = number_of_nodes - rank*number_of_local_nodes
    END IF
    ALLOCATE(local_nodes(number_of_local_nodes))
    fillrow: DO counter=1,number_of_local_nodes
       local_nodes(counter) = starting_node + (counter-1)
    END DO fillrow
    ending_node = local_nodes(number_of_local_nodes)
  END SUBROUTINE DivideUpWork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FillMatrix
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: num_of_edges
    INTEGER :: counter
    INTEGER, DIMENSION(:), ALLOCATABLE :: extra_scratch
    REAL :: temporary
    INTEGER :: extra_source_node, extra_destination_node

    CALL ConstructTripletList(triplet_list)

    !! First add the connection between each node and itself.
    fill_diagonal: DO counter=1, number_of_local_nodes
       temp_triplet%index_row = local_nodes(counter)
       temp_triplet%index_column = local_nodes(counter)
       temp_triplet%point_value = 1
       CALL AppendToTripletList(triplet_list,temp_triplet)
    END DO fill_diagonal

    !! Now connections between nearest neighbors.
    fill_neighbor: DO counter=1, number_of_local_nodes
       temp_triplet%index_row = local_nodes(counter)
       temp_triplet%point_value = 0.1
       IF (local_nodes(counter) .EQ. 1) THEN
          !! Right value
          temp_triplet%index_column = local_nodes(counter) + 1
          CALL AppendToTripletList(triplet_list,temp_triplet)
       ELSE IF (local_nodes(counter) .EQ. number_of_nodes) THEN
          !! Left value
          temp_triplet%index_column = local_nodes(counter) - 1
          CALL AppendToTripletList(triplet_list,temp_triplet)
       ELSE
          !! Left value
          temp_triplet%index_column = local_nodes(counter) - 1
          CALL AppendToTripletList(triplet_list,temp_triplet)
          !! Right value
          temp_triplet%index_column = local_nodes(counter) + 1
          CALL AppendToTripletList(triplet_list,temp_triplet)
       END IF
    END DO fill_neighbor

    !! Finally the random extra connections.
    ALLOCATE(extra_scratch(number_of_nodes))
    extra_scratch = 0
    counter = 1
    DO WHILE(counter .LE. extra_connections)
       CALL RANDOM_NUMBER(temporary)
       extra_source_node = CEILING(temporary*number_of_nodes)
       CALL RANDOM_NUMBER(temporary)
       extra_destination_node = CEILING(temporary*number_of_nodes)

       !! Check if we haven't already add this connection
       IF (extra_scratch(extra_source_node) .NE. 1 .AND. &
            & extra_scratch(extra_destination_node) .NE. 1 .AND. &
            & extra_source_node .NE. extra_destination_node .AND. &
            & extra_source_node .NE. extra_destination_node - 1 .AND. &
            & extra_source_node .NE. extra_destination_node + 1) THEN
          counter = counter + 1
          extra_scratch(extra_source_node) = 1
          extra_scratch(extra_destination_node) = 1

          IF (extra_source_node .GE. starting_node .AND. &
               & extra_source_node .LE. ending_node) THEN
             temp_triplet%index_row = extra_source_node
             temp_triplet%index_column = extra_destination_node
             temp_triplet%point_value = 0.1
             CALL AppendToTripletList(triplet_list, temp_triplet)
          END IF
          IF (extra_destination_node .GE. starting_node .AND. &
               & extra_destination_node .LE. ending_node) THEN
             temp_triplet%index_row = extra_destination_node
             temp_triplet%index_column = extra_source_node
             temp_triplet%point_value = 0.1
             CALL AppendToTripletList(triplet_list, temp_triplet)
          END IF
       END IF
    END DO

    !! Finally build the matrix
    CALL FillMatrixFromTripletList(NetworkMat, triplet_list)

    !! Cleanup
    DEALLOCATE(extra_scratch)
  END SUBROUTINE FillMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SolveMatrix
    TYPE(Matrix_ps) :: ResMat

    !! Compute attenuation*Identity - Matrix
    CALL ConstructEmptyMatrix(ResMat, number_of_nodes)
    CALL FillMatrixIdentity(ResMat)
    CALL IncrementMatrix(NetworkMat, ResMat, &
         & alpha_in=REAL(-1.0*attenuation,NTREAL))

    !! Invert
    CALL Invert(ResMat, ResultMat, solver_parameters)

    !! Cleanup
    CALL DestructMatrix(ResMat)
  END SUBROUTINE SolveMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM GraphTheory
