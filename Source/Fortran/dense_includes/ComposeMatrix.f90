  !! Local Data
  INTEGER, DIMENSION(block_rows + 1) :: row_offsets
  INTEGER, DIMENSION(block_columns + 1) :: column_offsets
  INTEGER :: out_rows, out_columns
  INTEGER :: II, JJ

  !! Determine the size of the big matrix
  out_columns = 0
  column_offsets(1) = 1
  out_columns = 0
  DO JJ = 1, block_columns
     column_offsets(JJ + 1) = column_offsets(JJ) + mat_array(1, JJ)%columns
     out_columns = out_columns + mat_array(1, JJ)%columns
  END DO
  row_offsets(1) = 1
  out_rows = 0
  DO II = 1, block_rows
     row_offsets(II + 1) = row_offsets(II) + mat_array(II, 1)%rows
     out_rows = out_rows + mat_array(II, 1)%rows
  END DO

  !! Allocate Memory
  CALL ConstructEmptyMatrix(out_matrix, out_columns, out_rows)

  !! Copy
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        out_matrix%DATA(row_offsets(II):row_offsets(II + 1) - 1, &
             & column_offsets(JJ):column_offsets(JJ + 1) - 1) = &
             & mat_array(II, JJ)%DATA
     END DO
  END DO
