  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  INTEGER :: rows, cols
  INTEGER II, JJ

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0_NTREAL
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0_NTREAL
  ELSE
     threshold = threshold_in
  END IF

  rows = GetMatrixBlockRows(matA)
  cols = GetMatrixBlockColumns(matA)

  !! Check that the target matrix is allocated.
  IF (.NOT. IsAllocated(matB)) THEN
     CALL ConstructEmptyMatrix(matB, rows, cols, IsBaseCase(matB))
  END IF
 
  IF (IsBaseCase(matA)) THEN
     DO II = 1, rows
       DO JJ = 1, cols
          CALL IncrementMatrix(matA%base_data(II,JJ), matB%base_data(II,JJ), &
               & alpha, threshold)
       END DO
     END DO
  ELSE 
     DO II = 1, rows
       DO JJ = 1, cols
          CALL IncrementMatrix(matA%h_data(II,JJ), matB%h_data(II,JJ), &
               & alpha, threshold)
       END DO
     END DO
  END IF