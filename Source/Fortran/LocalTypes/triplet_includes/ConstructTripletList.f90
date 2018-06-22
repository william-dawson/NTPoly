  CALL this%Destruct

  IF (PRESENT(size_in)) THEN
     this%CurrentSize  = size_in
  ELSE
     this%CurrentSize  = 0
  END IF

  ALLOCATE(this%data(this%CurrentSize))
