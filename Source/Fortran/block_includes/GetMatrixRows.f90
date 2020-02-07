  INTEGER :: rows
  INTEGER :: temp
  INTEGER II
  
  rows = GetMatrixBlockRows(this)

  IF (IsBaseCase(this)) THEN
     DO II = 1, rows
        total_rows = total_rows + this%base_data(II,1)%rows
     END DO
  ELSE
     DO II = 1, rows
     	temp = GetMatrixRows(this%h_data(II,1))
     	total_rows = total_rows + temp
     END DO
  END IF