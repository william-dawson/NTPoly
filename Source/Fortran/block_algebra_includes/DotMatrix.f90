  INTEGER :: rows, cols
  INTEGER II, JJ
  
  rows = GetMatrixBlockRows(matA)
  cols = GetMatrixBlockColumns(matA)

  product = 0.0_NTREAL
 
  IF (IsBaseCase(matA)) THEN
     DO II = 1, rows
        DO JJ = 1, cols
           CALL DotMatrix(matA%base_data(II,JJ), matB%base_data(II,JJ), temp)
           product = product + temp
        END DO
     END DO
  ELSE
     DO II = 1, rows
        DO JJ = 1, cols
           CALL DotMatrix(matA%h_data(II,JJ), matB%h_data(II,JJ), temp)
           product = product + temp
        END DO
     END DO
  END IF