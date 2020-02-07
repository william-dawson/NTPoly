  INTEGER :: rows, cols
  INTEGER :: II, JJ

  IF (IsBaseCase(this)) THEN
     CALL ComposeMatrix(this%base_data, out_matrix)
  ELSE 
     rows = GetMatrixBlockRows(this)
     cols = GetMatrixBlockColumns(this)
     CALL ConstructEmptyMatrix(base_mat, rows, cols, base_case=.TRUE.)
     DO II = 1, rows
       DO JJ = 1, cols
         CALL ComposeMatrix(this%h_data(II,JJ), base_mat%base_data(II,JJ))
       END DO
     END DO
     !! Then compose the base case version.
     CALL ComposeMatrix(base_mat, out_matrix)
     CALL DestructMatrix(base_mat)
  END IF