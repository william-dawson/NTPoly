  INTEGER, INTENT(IN) :: rank
  TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
  INTEGER, INTENT(IN) :: send_tag
  INTEGER :: comm
  !! Local Data
  INTEGER :: ierr

  !! Send the data
  CALL MPI_Isend(inmat%inner_index, SIZE(inmat%values), MPI_INT, &
       & rank, send_tag, comm, helper%inner_request, ierr)
  CALL MPI_Isend(inmat%outer_index, inmat%columns+1, MPI_INT, &
       & rank, send_tag, comm, helper%outer_request, ierr)
