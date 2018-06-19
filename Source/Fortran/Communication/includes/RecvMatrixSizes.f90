  INTEGER, INTENT(IN) :: rank
  INTEGER, INTENT(INOUT) :: comm
  TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
  INTEGER, INTENT(IN) :: recv_tag
  !! Local Data
  INTEGER :: ierr

  CALL MPI_Irecv(helper%matrix_data, 3, MPI_INTEGER, rank, recv_tag, comm, &
       & helper%size_request, ierr)

 helper%recv_stage = 0
