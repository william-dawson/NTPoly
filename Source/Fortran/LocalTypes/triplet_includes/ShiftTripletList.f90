  !! Local Variables
  INTEGER :: II

  DO II = 1, this%CurrentSize
     this%data(II)%index_row = this%data(II)%index_row + row_shift
     this%data(II)%index_column = this%data(II)%index_column + column_shift
  END DO
