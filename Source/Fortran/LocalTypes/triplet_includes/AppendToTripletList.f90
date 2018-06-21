  old_size = SIZE(this%data)

  !! First, check if we need to allocate more memory
  IF (this%CurrentSize+1 .GT. old_size) THEN
     IF (old_size .EQ. 0) THEN
        new_size = 1
     ELSE IF (old_size .EQ. 1) THEN
        new_size = 2
     ELSE
        new_size = INT(old_size*1.5)
     END IF
     CALL this%Resize(new_size)
  END IF

  !! Append
  this%CurrentSize = this%CurrentSize+1

  this%data(this%CurrentSize) = triplet_value
