  CALL this%Destruct

  this%rows = rows
  this%columns = columns

  !! Allocate
  ALLOCATE(this%outer_index(this%columns+1))
  this%outer_index = 0
  ALLOCATE(this%inner_index(triplet_list%CurrentSize))
  ALLOCATE(this%values(triplet_list%CurrentSize))

  !! Insert Values
  outer_array_ptr = 1
  DO values_counter = 1, triplet_list%CurrentSize
     !! Moving on to the next column?
     DO WHILE(.NOT. triplet_list%data(values_counter)%index_column .EQ. &
          & outer_array_ptr)
        outer_array_ptr = outer_array_ptr + 1
        this%outer_index(outer_array_ptr+1) = this%outer_index(outer_array_ptr)
     END DO
     this%outer_index(outer_array_ptr+1)=this%outer_index(outer_array_ptr+1)+1
     !! Insert inner index and value
     this%inner_index(values_counter) = &
          & triplet_list%data(values_counter)%index_row
     this%values(values_counter) = &
          & triplet_list%data(values_counter)%point_value
  END DO

  !! Fill In The Rest Of The Outer Values
  DO outer_array_ptr = outer_array_ptr+2, this%columns+1
     this%outer_index(outer_array_ptr) = this%outer_index(outer_array_ptr-1)
  END DO
