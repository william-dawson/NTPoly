  !! Local Variables
  INTEGER :: pind
  INTEGER :: II
  INTEGER :: swap
  INTEGER :: ierr

  !! Search for the maximum pivot
  max_diag = [0, 0]
  DO II = start_index, AMat%actual_matrix_dimension
     pind = pivot_vector(II)
     IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
        temp_diag = diag(pind - AMat%start_column + 1)
        IF (temp_diag .GT. max_diag(1)) THEN
           max_diag(1) = temp_diag
           max_diag(2) = II
        END IF
     END IF
  END DO
  CALL MPI_Allreduce(MPI_IN_PLACE, max_diag, 1, MPI_2DOUBLE_PRECISION, &
       & MPI_MAXLOC, process_grid%row_comm, ierr)
  index = INT(max_diag(2))
  value = max_diag(1)

  swap = pivot_vector(index)
  pivot_vector(index) = pivot_vector(start_index)
  pivot_vector(start_index) = swap

  index = swap

  !! Determine local pivots
  num_local_pivots = 0
  DO II = start_index + 1, AMat%actual_matrix_dimension
     pind = pivot_vector(II)
     IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
        num_local_pivots = num_local_pivots + 1
        local_pivots(num_local_pivots) = pind - AMat%start_column + 1
     END IF
  END DO
