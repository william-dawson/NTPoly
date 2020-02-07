  IF (IsBaseCase(this)) THEN
   	 rows = SIZE(this%base_data, DIM=1)
  ELSE
     rows = SIZE(this%h_data, DIM=1)
  END IF