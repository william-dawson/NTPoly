  INTEGER :: rows, cols
  INTEGER II, JJ

  rows = GetMatrixBlockRows(matA)
  cols = GetMatrixBlockColumns(matA)

  IF (IsBaseCase(matA)) THEN
     CALL ConstructEmptyMatrix(matB, rows, cols, .TRUE.)
     DO II = 1, rows
       DO JJ = 1, cols
          CALL CopyMatrix(matA%base_data(II,JJ), matB%base_data(II,JJ))
       END DO
     END DO
  ELSE 
     CALL ConstructEmptyMatrix(matB, rows, cols, .FALSE.)
     DO II = 1, rows
       DO JJ = 1, cols
          CALL CopyMatrix(matA%h_data(II,JJ), matB%h_data(II,JJ))
       END DO
     END DO
  END IF