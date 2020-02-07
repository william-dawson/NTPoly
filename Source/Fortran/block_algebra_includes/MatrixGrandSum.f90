  INTEGER :: rows, cols
  INTEGER II, JJ
  
  rows = GetMatrixBlockRows(this)
  cols = GetMatrixBlockColumns(this)

  sum_value = 0.0_NTREAL

  !! Hierarchical case. 
  IF (IsBaseCase(this)) THEN
     DO II = 1, rows
        DO JJ = 1, cols
           CALL MatrixGrandSum(this%base_data(II,JJ), temp)
           sum_value = sum_value + temp
        END DO
     END DO
  ELSE
     DO II = 1, rows
        DO JJ = 1, cols
           CALL MatrixGrandSum(this%h_data(II,JJ), temp)
           sum_value = sum_value + temp
        END DO
     END DO
  END IF