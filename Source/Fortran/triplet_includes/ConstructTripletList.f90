  !! Local data
  INTEGER :: size

  IF (PRESENT(size_in)) THEN
     size = size_in
  ELSE
     size = 0
  END IF

  IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
  this%CurrentSize = size

  ALLOCATE(this%DATA(size))
