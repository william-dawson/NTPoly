  !! Local Data
  INTEGER, DIMENSION(block_rows) :: block_size_row
  INTEGER, DIMENSION(block_columns) :: block_size_column
  INTEGER, DIMENSION(block_rows+1) :: row_offsets
  INTEGER, DIMENSION(block_columns+1) :: column_offsets
  !! Temporary Variables
  INTEGER :: divisor_row, divisor_column
  INTEGER :: II, JJ

  !! Calculate the split sizes
  IF (PRESENT(block_size_row_in)) THEN
     block_size_row = block_size_row_in
  ELSE
     divisor_row = this%rows/block_rows
     block_size_row = divisor_row
     block_size_row(block_rows) = this%rows - divisor_row*(block_rows-1)
  END IF
  IF (PRESENT(block_size_column_in)) THEN
     block_size_column = block_size_column_in
  ELSE
     divisor_column = this%columns/block_columns
     block_size_column = divisor_column
     block_size_column(block_columns) = this%columns - &
          & divisor_column*(block_columns-1)
  END IF

  !! Copy the block offsets
  row_offsets(1) = 1
  DO II = 1, block_rows
     row_offsets(II+1) = row_offsets(II) + block_size_row(II)
  END DO
  column_offsets(1) = 1
  DO JJ = 1, block_columns
     column_offsets(JJ+1) = column_offsets(JJ) + block_size_column(JJ)
  END DO

  !! Copy
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL ConstructEmptyMatrix(split_array(II,JJ), block_size_column(JJ), &
             & block_size_row(II))
        split_array(II,JJ)%DATA = &
             & this%DATA(row_offsets(II):row_offsets(II+1)-1, &
             & column_offsets(JJ):column_offsets(JJ+1)-1)
     END DO
  END DO
