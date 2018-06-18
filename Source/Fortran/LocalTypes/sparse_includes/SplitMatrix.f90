  !! Local Data
  INTEGER, DIMENSION(block_rows) :: block_size_row
  INTEGER, DIMENSION(block_columns) :: block_size_column
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

  !! First split by columns which is easy with the CSR format
  CALL StartTimer("Split Column")
  CALL SplitMatrixColumns(this, block_columns, block_size_column, &
       & column_split)
  CALL StopTimer("Split Column")

  !! Now Split By Rows
  CALL StartTimer("Split Row")
  DO JJ = 1, block_columns
     CALL TransposeMatrix(column_split(JJ), Temp)
     CALL SplitMatrixColumns(Temp, block_rows, block_size_row, &
          & row_split)
     !! Copy into output array
     DO II = 1, block_rows
        CALL TransposeMatrix(row_split(II), split_array(II,JJ))
     END DO
  END DO
  CALL StopTimer("Split Row")

  !! Cleanup
  CALL DestructMatrix(Temp)
  DO II = 1, block_rows
     CALL DestructMatrix(row_split(II))
  END DO
  DO II = 1, block_columns
     CALL DestructMatrix(column_split(II))
  END DO
