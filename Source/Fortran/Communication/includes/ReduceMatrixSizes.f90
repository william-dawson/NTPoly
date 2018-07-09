  !! Local Data
  INTEGER :: grid_error

  CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

  !! Gather Information About Other Processes
  ALLOCATE(helper%values_per_process(helper%comm_size))
  CALL MPI_IAllGather(SIZE(matrix%values),1,MPI_INT,&
       & helper%values_per_process,1,MPI_INT,communicator, &
       & helper%size_request, grid_error)
