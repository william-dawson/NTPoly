  !! Local data
  INTEGER :: size

  IF (PRESENT(size_in)) THEN
     size = size_in
  ELSE
     size = 0
  END IF

  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  this%CurrentSize = size

  ALLOCATE(this%data(size))
