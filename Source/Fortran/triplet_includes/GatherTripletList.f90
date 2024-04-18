!! Local Data - Send/Recv Buffers
INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col

!! Sizes help
INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts
INTEGER, DIMENSION(:), ALLOCATABLE :: displ
INTEGER :: gather_size

!! Temporary variables
INTEGER :: num_processes, II, ierr

!! Figure out the comm size
CALL MPI_COMM_SIZE(comm, num_processes, ierr)

!! Get the count 
ALLOCATE(recvcounts(num_processes))
CALL MPI_Allgather(triplet_in%CurrentSize, 1, MPI_INTEGER, recvcounts, &
     & 1, MPI_INTEGER, comm, ierr)

!! Get the displacements
gather_size = SUM(recvcounts)
ALLOCATE(displ(num_processes))
displ(1) = 0
DO II = 2, num_processes
  displ(II) = displ(II - 1) + recvcounts(II - 1)
END DO

!! Prepare the send buffers
ALLOCATE(send_buffer_row(triplet_in%CurrentSize))
ALLOCATE(send_buffer_col(triplet_in%CurrentSize))
ALLOCATE(send_buffer_val(triplet_in%CurrentSize))
DO II = 1, triplet_in%CurrentSize
   CALL GetTripletAt(triplet_in, II, temp_triplet)
   send_buffer_row(II) = temp_triplet%index_row
   send_buffer_col(II) = temp_triplet%index_column
   send_buffer_val(II) = temp_triplet%point_value
END DO

!! Gather Call
ALLOCATE(recv_buffer_row(gather_size))
ALLOCATE(recv_buffer_col(gather_size))
ALLOCATE(recv_buffer_val(gather_size))
CALL MPI_Allgatherv(send_buffer_row, triplet_in%CurrentSize, MPI_INTEGER, &
     & recv_buffer_row, recvcounts, displ, MPI_INTEGER, comm, ierr)
CALL MPI_Allgatherv(send_buffer_col, triplet_in%CurrentSize, MPI_INTEGER, &
     & recv_buffer_col, recvcounts, displ, MPI_INTEGER, comm, ierr)
CALL MPI_Allgatherv(send_buffer_val, triplet_in%CurrentSize, MPIDATATYPE, &
     & recv_buffer_val, recvcounts, displ, MPIDATATYPE, comm, ierr)

!! Unpack
CALL ConstructTripletList(gathered_out, gather_size)
DO II = 1, gather_size
   gathered_out%DATA(II)%index_row = recv_buffer_row(II)
   gathered_out%DATA(II)%index_column = recv_buffer_col(II)
   gathered_out%DATA(II)%point_value = recv_buffer_val(II)
END DO

!! Cleanup
DEALLOCATE(recvcounts)
DEALLOCATE(displ)
DEALLOCATE(send_buffer_row)
DEALLOCATE(send_buffer_col)
DEALLOCATE(send_buffer_val)
