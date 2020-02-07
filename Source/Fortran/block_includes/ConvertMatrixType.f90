  INTEGER :: rows, cols
  INTEGER :: II, JJ

  rows = GetMatrixBlockRows(inmat)
  cols = GetMatrixBlockColumns(inmat)

  IF (IsBaseCase(inmat)) THEN
     DO II = 1, rows
       DO JJ = 1, cols
       	 CALL ConstructMatrix(outmat, rows, cols, .TRUE.)
         CALL ConvertMatrixType(inmat%base_data(II,JJ), &
         	  & outmat%base_data(II,JJ))
       END DO
     END DO
  ELSE
     DO II = 1, rows
       DO JJ = 1, cols
       	 CALL ConstructMatrix(outmat, rows, cols, .FALSE.)
         CALL ConvertMatrixType(inmat%h_data(II,JJ), outmat%h_data(II,JJ))
       END DO
     END DO
  END IF