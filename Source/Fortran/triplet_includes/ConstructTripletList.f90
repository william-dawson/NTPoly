  !! Local data
  INTEGER :: size

  IF (PRESENT(size_in)) THEN
     size = size_in
  ELSE
     size = 0
  END IF

  CALL DestructTripletList(this)

  this%CurrentSize = size

  ALLOCATE(this%DATA(size))
