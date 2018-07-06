  !! Local data
  INTEGER :: new_size
  INTEGER :: idata

  DO idata = 1, this%CurrentSize
     IF ((this%DATA(idata)%index_row == triplet_value%index_row) .AND. &
          & (this%DATA(idata)%index_column == triplet_value%index_column)) THEN
        this%DATA(idata)%point_value = this%DATA(idata)%point_value + &
             & triplet_value%point_value
        RETURN
     END IF
  END DO

  !! First, check if we need to allocate more memory
  IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
     IF (SIZE(this%data) .EQ. 0) THEN
        new_size = 1
     ELSE IF (SIZE(this%data) .EQ. 1) THEN
        new_size = 2
     ELSE
        new_size = INT(SIZE(this%data)*1.5)
     END IF
     CALL ResizeTripletList(this,new_size)
  END IF

  !! Append
  this%CurrentSize = this%CurrentSize+1
  this%DATA(this%CurrentSize) = triplet_value
