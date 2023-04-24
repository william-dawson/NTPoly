  !! Helper Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: elements_per_inner
  INTEGER :: total_counter

  CALL ConstructEmptyMatrix(dense_matrix, sparse_matrix%rows, &
       & sparse_matrix%columns)

  !! Loop over elements.
  dense_matrix%DATA = 0
  total_counter = 1
  DO outer_counter = 1, sparse_matrix%columns
     elements_per_inner = sparse_matrix%outer_index(outer_counter + 1) - &
          & sparse_matrix%outer_index(outer_counter)
     temporary%index_column = outer_counter
     DO inner_counter = 1, elements_per_inner
        temporary%index_row = sparse_matrix%inner_index(total_counter)
        temporary%point_value = sparse_matrix%values(total_counter)
        dense_matrix%DATA(temporary%index_row, temporary%index_column) = &
             & temporary%point_value
        total_counter = total_counter + 1
     END DO
  END DO
