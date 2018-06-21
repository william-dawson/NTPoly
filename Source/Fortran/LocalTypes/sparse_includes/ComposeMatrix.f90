  INTEGER :: II, JJ

  !! First transpose the matrices
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL mat_t(II,JJ)%Transpose(mat_array(II,JJ))
     END DO
  END DO

  !! Next merge the columns
  DO JJ = 1, block_columns
     CALL ComposeMatrixColumns(mat_t(:,JJ), Temp)
     CALL merged_columns(JJ)%Transpose(Temp)
  END DO

  !! Final Merge
  CALL ComposeMatrixColumns(merged_columns, out_matrix)

  !! Cleanup
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL mat_t(II,JJ)%Destruct
     END DO
  END DO
  DO JJ = 1, block_columns
     CALL merged_columns(JJ)%Destruct
  END DO
