  !! Loop
  DO counter = 1, triplet_list%CurrentSize
     triplet_list%data(counter)%index_row = &
          triplet_list%data(counter)%index_row + row_shift
     triplet_list%data(counter)%index_column = &
          triplet_list%data(counter)%index_column + column_shift
  END DO
