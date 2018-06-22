  this%rows = rows
  this%columns = columns

  ALLOCATE(this%data(rows,columns))

  IF (PRESENT(zero_in)) THEN
    IF (zero_in) THEN
      this%data = 0
    END IF
  END IF
