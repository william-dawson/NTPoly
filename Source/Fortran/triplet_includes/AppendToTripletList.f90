  !! Local data
  INTEGER :: new_size

  !! First, check if we need to allocate more memory
  IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
     IF (SIZE(this%DATA) .EQ. 0) THEN
        new_size = 1
     ELSE IF (SIZE(this%DATA) .EQ. 1) THEN
        new_size = 2
     ELSE
        new_size = INT(SIZE(this%DATA)*1.5)
     END IF
     CALL ResizeTripletList(this,new_size)
  END IF

  !! Append
  this%CurrentSize = this%CurrentSize+1
  this%DATA(this%CurrentSize) = triplet_value
