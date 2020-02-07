  INTEGER :: cols
  INTEGER :: temp
  INTEGER JJ
  
  cols = GetMatrixBlockColumns(this)

  IF (IsBaseCase(this)) THEN
     DO JJ = 1, cols
        total_columns = total_columns + this%base_data(1,JJ)%columns
     END DO
  ELSE
     DO JJ = 1, cols
     	temp = GetMatrixColumns(this%h_data(1,JJ))
     	total_columns = total_columns + temp
     END DO
  END IF