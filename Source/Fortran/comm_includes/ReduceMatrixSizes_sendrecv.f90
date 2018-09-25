  !! Local Data
  INTEGER :: grid_error
  INTEGER :: sum_outer_indices
  INTEGER :: II

  CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

  !! Build Storage
  CALL DestructMatrix(gathered_matrix)
  sum_outer_indices = (matrix%columns+1)*helper%comm_size
  ALLOCATE(gathered_matrix%outer_index(sum_outer_indices+1))

 !! Send/Recv Outer Index
  DO II = 1, helper%comm_size
     !! Send/Recv Outer Index
     CALL MPI_ISend(matrix%outer_index, SIZE(matrix%outer_index), MPI_INT, &
          & II-1, 3, communicator, helper%outer_send_request_list(II), &
          & grid_error)
     istart = (matrix%columns+1)*(II-1)+1
     isize = matrix%columns + 1
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%outer_index(istart:iend), isize, MPI_INT, &
          & II-1, 3, communicator, &
          & helper%outer_recv_request_list(II), grid_error)
  END DO
