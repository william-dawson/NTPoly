  !! Local Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: columns, rows

  columns = dense_matrix%columns
  rows = dense_matrix%rows

  IF (PRESENT(threshold_in)) THEN
     CALL temporary_list%Init
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           IF (ABS(temporary%point_value) .GT. threshold_in) THEN
              temporary%index_row = inner_counter
              CALL temporary_list%Append(temporary)
           END IF
        END DO
     END DO
  ELSE
     CALL temporary_list%Init(rows*columns)
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           temporary%index_row = inner_counter
           temporary_list%data(inner_counter+rows*(outer_counter-1)) = &
                & temporary
        END DO
     END DO
  END IF

  CALL sparse_matrix%InitFromTripletList(temporary_list, rows, columns)
