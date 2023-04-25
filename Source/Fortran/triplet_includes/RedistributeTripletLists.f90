  !! Local Data - Offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_offsets
  !! Local Data - Send/Recv Buffers
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
  !! ETC
  INTEGER :: num_processes
  INTEGER :: counter, inner_counter, insert_pt
  INTEGER :: mpi_error

  !! Allocate Size Buffers
  CALL MPI_COMM_SIZE(comm, num_processes, mpi_error)
  ALLOCATE(send_per_process(num_processes))
  ALLOCATE(send_offsets(num_processes))
  ALLOCATE(recv_per_process(num_processes))
  ALLOCATE(recv_offsets(num_processes))

  !! Figure Out How Much Data Gets Sent
  DO counter = 1, num_processes
     send_per_process(counter) = triplet_lists(counter)%CurrentSize
  END DO
  send_offsets(1) = 0
  DO counter = 2, num_processes
     send_offsets(counter) = send_offsets(counter - 1) + &
          & send_per_process(counter - 1)
  END DO

  !! Figure Out How Much Data Gets Received
  CALL MPI_ALLTOALL(send_per_process, 1, MPINTINTEGER, recv_per_process, 1, &
       & MPINTINTEGER, comm, mpi_error)
  recv_offsets(1) = 0
  DO counter = 2, num_processes
     recv_offsets(counter) = recv_offsets(counter - 1) + &
          & recv_per_process(counter - 1)
  END DO

  !! Allocate And Fill Send Buffers
  ALLOCATE(send_buffer_row(SUM(send_per_process)))
  ALLOCATE(send_buffer_col(SUM(send_per_process)))
  ALLOCATE(send_buffer_val(SUM(send_per_process)))
  ALLOCATE(recv_buffer_row(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_col(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_val(SUM(recv_per_process)))

  !! Fill Send Buffer
  insert_pt = 1
  DO counter = 1, num_processes
     DO inner_counter = 1, triplet_lists(counter)%CurrentSize
        CALL GetTripletAt(triplet_lists(counter), inner_counter, temp_triplet)
        send_buffer_row(insert_pt) = temp_triplet%index_row
        send_buffer_col(insert_pt) = temp_triplet%index_column
        send_buffer_val(insert_pt) = temp_triplet%point_value
        insert_pt = insert_pt + 1
     END DO
  END DO

  !! Do Actual Send
  CALL MPI_Alltoallv(send_buffer_col, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_col, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_row, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_row, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_val, send_per_process, send_offsets, &
       & MPIDATATYPE, recv_buffer_val, recv_per_process, recv_offsets, &
       & MPIDATATYPE, comm, mpi_error)

  !! Unpack Into The Output Triplet List
  CALL ConstructTripletList(local_data_out, size_in=SUM(recv_per_process))
  DO counter = 1, SUM(recv_per_process)
     local_data_out%DATA(counter)%index_column = recv_buffer_col(counter)
     local_data_out%DATA(counter)%index_row = recv_buffer_row(counter)
     local_data_out%DATA(counter)%point_value = recv_buffer_val(counter)
  END DO

  !! Cleanup
  DEALLOCATE(send_per_process)
  DEALLOCATE(send_offsets)
  DEALLOCATE(recv_per_process)
  DEALLOCATE(recv_offsets)
  DEALLOCATE(send_buffer_row)
  DEALLOCATE(send_buffer_col)
  DEALLOCATE(send_buffer_val)
  DEALLOCATE(recv_buffer_row)
  DEALLOCATE(recv_buffer_col)
  DEALLOCATE(recv_buffer_val)
