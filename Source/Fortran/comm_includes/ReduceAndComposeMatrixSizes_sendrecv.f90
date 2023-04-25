  !! Local Data
  INTEGER :: II
  INTEGER :: istart, isize, iend
  INTEGER :: ierr

  CALL MPI_Comm_size(comm, helper%comm_size, ierr)

  !! Build Storage
  CALL ConstructEmptyMatrix(gathered_matrix, &
       & matrix%rows, matrix%columns * helper%comm_size)
  gathered_matrix%outer_index(1) = 0

  ALLOCATE(helper%outer_send_request_list(helper%comm_size))
  ALLOCATE(helper%outer_recv_request_list(helper%comm_size))

  !! Send/Recv Outer Index
  DO II = 1, helper%comm_size
     !! Send/Recv Outer Index
     CALL MPI_ISend(matrix%outer_index(2:), matrix%columns, MPINTINTEGER, &
          & II-1, 3, comm, helper%outer_send_request_list(II), ierr)
     istart = (matrix%columns)*(II - 1) + 2
     isize = matrix%columns
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%outer_index(istart:iend), isize, &
          & MPINTINTEGER, II - 1, 3, comm, &
          & helper%outer_recv_request_list(II), ierr)
  END DO
