  
  INTEGER :: col, II

  DO II = 1, tlist%CurrentSize
     col = tlist%DATA(II)%index_column
     val = tlist%DATA(II)%point_value
     mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1)) = &
          val * mat%values(mat%outer_index(col) + 1:mat%outer_index(col + 1))
  END DO
