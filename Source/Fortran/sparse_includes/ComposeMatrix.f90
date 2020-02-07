  INTEGER :: block_columns, block_rows
  INTEGER :: II, JJ

  block_rows = SIZE(mat_array, DIM=1)
  block_columns = SIZE(mat_array, DIM=2)

  ALLOCATE(merged_columns(block_columns))
  ALLOCATE(mat_t(block_rows,block_columns))

  !! First transpose the matrices
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL TransposeMatrix(mat_array(II,JJ), mat_t(II,JJ))
     END DO
  END DO

  !! Next merge the columns
  DO JJ = 1, block_columns
     CALL ComposeMatrixColumns(mat_t(:,JJ), Temp)
     CALL TransposeMatrix(Temp, merged_columns(JJ))
  END DO

  !! Final Merge
  CALL ComposeMatrixColumns(merged_columns, out_matrix)

  !! Cleanup
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL DestructMatrix(mat_t(II,JJ))
     END DO
  END DO
  DO JJ = 1, block_columns
     CALL DestructMatrix(merged_columns(JJ))
  END DO

  DEALLOCATE(merged_columns)
  DEALLOCATE(mat_t)