  !! Local Variables
  INTEGER :: fill_counter
  INTEGER :: II, JJ, local_JJ, local_II
  INTEGER, DIMENSION(process_grid%num_process_rows) :: diags_per_proc
  INTEGER, DIMENSION(process_grid%num_process_rows) :: diag_displ
  INTEGER :: ierr

  diag = 0
  fill_counter = 0
  DO JJ = AMat%start_column, AMat%end_column - 1
     IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
        local_JJ = JJ - AMat%start_column + 1
        local_II = JJ - AMat%start_row + 1
        diag(local_JJ) = dense_a%DATA(local_II, local_JJ)
        fill_counter = fill_counter + 1
     END IF
  END DO
  diags_per_proc(process_grid%my_row + 1) = fill_counter

  !! Duplicate the diagonal entries along the process column (across rows)
  CALL MPI_Allgather(MPI_IN_PLACE, 1, MPINTINTEGER, diags_per_proc, 1, &
       & MPINTINTEGER, process_grid%column_comm, ierr)
  diag_displ(1) = 0
  DO II = 2, process_grid%num_process_rows
     diag_displ(II) = diag_displ(II - 1) + diags_per_proc(II - 1)
  END DO
