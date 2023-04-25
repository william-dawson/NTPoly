  !! Local Data
  INTEGER :: II, outer_idx

  CALL DestructMatrix(this)

  this%rows = rows
  this%columns = columns

  !! Allocate
  ALLOCATE(this%outer_index(this%columns + 1))
  this%outer_index = 0
  ALLOCATE(this%inner_index(triplet_list%CurrentSize))
  ALLOCATE(this%values(triplet_list%CurrentSize))

  !! Insert Values
  outer_idx = 1
  DO II = 1, triplet_list%CurrentSize
     !! Moving on to the next column?
     DO WHILE(.NOT. triplet_list%DATA(II)%index_column .EQ. outer_idx)
        outer_idx = outer_idx + 1
        this%outer_index(outer_idx + 1) = this%outer_index(outer_idx)
     END DO
     this%outer_index(outer_idx + 1)=this%outer_index(outer_idx + 1) + 1
     !! Insert inner index and value
     this%inner_index(II) = triplet_list%DATA(II)%index_row
     this%values(II) = triplet_list%DATA(II)%point_value
  END DO

  !! Fill In The Rest Of The Outer Values
  DO outer_idx = outer_idx + 2, this%columns + 1
     this%outer_index(outer_idx) = this%outer_index(outer_idx - 1)
  END DO
