  INTEGER :: rows, cols
  INTEGER II, JJ
  
  rows = GetMatrixBlockRows(matA)
  cols = GetMatrixBlockColumns(matA)
 
  IF (IsBaseCase(matA)) THEN
     CALL ConstructEmptyMatrix(matC, rows, cols, .TRUE.)
     DO II = 1, rows
        DO JJ = 1, cols
           CALL PairwiseMultiplyMatrix(matA%base_data(II,JJ), &
           	    & matB%base_data(II,JJ), matC%base_data(II,JJ))
        END DO
     END DO
  ELSE 
     CALL ConstructEmptyMatrix(matC, rows, cols, .FALSE.)
     DO II = 1, rows
        DO JJ = 1, cols
           CALL PairwiseMultiplyMatrix(matA%h_data(II,JJ), &
           	    & matB%h_data(II,JJ), matC%h_data(II,JJ))
        END DO
     END DO
  END IF