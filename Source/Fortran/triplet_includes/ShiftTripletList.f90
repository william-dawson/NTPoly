  !! Loop
  DO counter = 1, triplet_list%CurrentSize
     triplet_list%DATA(counter)%index_row = &
          triplet_list%DATA(counter)%index_row + row_shift
     triplet_list%DATA(counter)%index_column = &
          triplet_list%DATA(counter)%index_column + column_shift
  END DO
