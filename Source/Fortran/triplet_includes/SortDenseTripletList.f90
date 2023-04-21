  !! Local Data
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: dirty_buffer
  INTEGER :: list_length
  INTEGER :: row, col, ind
  INTEGER :: II, JJ

  !! Setup Memory
  ALLOCATE(value_buffer(matrix_rows, matrix_columns))
  ALLOCATE(dirty_buffer(matrix_rows, matrix_columns))
  value_buffer = 0
  dirty_buffer = 0
  list_length = input_list%CurrentSize
  CALL ConstructTripletList(sorted_list, list_length)

  !! Unpack
  DO II = 1, list_length
     row = input_list%DATA(II)%index_row
     col = input_list%DATA(II)%index_column
     value_buffer(row,col) = input_list%DATA(II)%point_value
     dirty_buffer(row,col) = 1
  END DO

  !! Repack
  ind = 1
  DO JJ = 1, matrix_columns
     DO II = 1, matrix_rows
        IF (dirty_buffer(II,JJ) .EQ. 1) THEN
           sorted_list%DATA(ind)%index_row = II
           sorted_list%DATA(ind)%index_column = JJ
           sorted_list%DATA(ind)%point_value = value_buffer(II, JJ)
           ind = ind + 1
        END IF
     END DO
  END DO

  !! Cleanup
  DEALLOCATE(value_buffer)
  DEALLOCATE(dirty_buffer)
