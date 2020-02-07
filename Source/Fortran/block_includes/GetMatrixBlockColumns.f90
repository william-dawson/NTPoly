  IF (IsBaseCase(this)) THEN
   	 columns = SIZE(this%base_data, DIM=2)
  ELSE
     columns = SIZE(this%h_data, DIM=2)
  END IF