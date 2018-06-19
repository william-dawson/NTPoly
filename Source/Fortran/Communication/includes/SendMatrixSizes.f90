  INTEGER, INTENT(IN) :: rank
  INTEGER, INTENT(INOUT) :: comm
  TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
  INTEGER, INTENT(IN) :: send_tag
  !! Local Data
  INTEGER :: ierr

  helper%matrix_data(1) = inmat%rows*inmat%columns
  helper%matrix_data(2) = inmat%rows
  helper%matrix_data(3) = inmat%columns

  CALL MPI_Isend(helper%matrix_data, 3, MPI_INTEGER, rank, send_tag, comm, &
       & helper%size_request, ierr)
