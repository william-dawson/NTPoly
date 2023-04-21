  !! Loop
  DO II = 1, triplet_list%CurrentSize
     triplet_list%DATA(II)%index_row = &
          triplet_list%DATA(II)%index_row + row_shift
     triplet_list%DATA(II)%index_column = &
          triplet_list%DATA(II)%index_column + column_shift
  END DO
