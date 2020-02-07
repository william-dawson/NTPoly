  INTEGER :: rows, cols
  INTEGER :: II, JJ

  rows = GetMatrixBlockRows(this)
  cols = GetMatrixBlockColumns(this)

  !! Hierarchical case. 
  IF (ALLOCATED(this%h_data)) THEN
     DO II = 1, rows
        DO JJ = 1, cols
          CALL ScaleMatrix(this%h_data(II,JJ), constant)
        END DO
     END DO
  ELSE !! Back case
     DO II = 1, rows
        DO JJ = 1, cols
           CALL ScaleMatrix(this%base_data(II,JJ), constant)
        END DO
     END DO
  END IF