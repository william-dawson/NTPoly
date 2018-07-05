  this%rows = rows
  this%columns = columns
  ALLOCATE(this%outer_index(this%columns+1))
  this%outer_index = 0

  IF (PRESENT(zero_in)) THEN
     IF (zero_in) THEN
        ALLOCATE(this%inner_index(0))
        ALLOCATE(this%values(0))
     END IF
  END IF
