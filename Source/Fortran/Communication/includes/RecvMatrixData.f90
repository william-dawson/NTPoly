
  INTEGER, INTENT(IN) :: rank
  TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
  INTEGER, INTENT(IN) :: recv_tag
  INTEGER :: comm
  !! Local Data
  INTEGER :: ierr

  !! Build Storage
  CALL ConstructEmptyMatrix(outmat, helper%matrix_data(3), &
       & helper%matrix_data(2))
  ! outmat =  Matrix_lsr(helper%matrix_data(3), helper%matrix_data(2))
  ALLOCATE(outmat%values(helper%matrix_data(1)))
  ALLOCATE(outmat%inner_index(helper%matrix_data(1)))
  outmat%outer_index(1) = 0

  !! Receive the data
  CALL MPI_Irecv(outmat%inner_index, helper%matrix_data(1), MPI_INT, &
       & rank, recv_tag, comm, helper%inner_request, ierr)
  CALL MPI_Irecv(outmat%outer_index, helper%matrix_data(3)+1, MPI_INT, &
       & rank, recv_tag, comm, helper%outer_request, ierr)
