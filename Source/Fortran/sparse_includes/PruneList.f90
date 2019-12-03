  !! Local data
  INTEGER :: row_counter_c, column_counter_c, hash_counter
  INTEGER :: working_column
  INTEGER :: temp_values_per_hash
  INTEGER :: pruned_counter

  pruned_counter = 1
  DO row_counter_c = 1, mat_c_rows
     DO column_counter_c = 1, (mat_c_columns-1)/memorypool%hash_size+1
        !! Sort the elements in a hash
        temp_values_per_hash = memorypool%inserted_per_bucket(&
             & column_counter_c,row_counter_c)
        memorypool%inserted_per_bucket(column_counter_c,row_counter_c) = 0
        !! Copy them
        DO hash_counter=1,temp_values_per_hash
           working_column = memorypool%hash_index(hash_counter+ &
                & (column_counter_c-1)*memorypool%hash_size, row_counter_c)
           working_value = memorypool%value_array(working_column,row_counter_c)
           memorypool%value_array(working_column,row_counter_c) = 0
           memorypool%dirty_array(working_column,row_counter_c) = .FALSE.
           IF (ABS(alpha*working_value) .GT. threshold) THEN
              memorypool%pruned_list(pruned_counter)%point_value = &
                   & alpha*working_value
              memorypool%pruned_list(pruned_counter)%index_column = &
                   & working_column
              memorypool%pruned_list(pruned_counter)%index_row = &
                   & row_counter_c
              pruned_counter = pruned_counter + 1
           END IF
        END DO
     END DO
  END DO
  CALL ConstructTripletList(unsorted_pruned_list, pruned_counter-1)
  unsorted_pruned_list%data = memorypool%pruned_list(1:pruned_counter-1)
  CALL SortTripletList(unsorted_pruned_list, mat_c_columns, mat_c_rows, &
       & sorted_pruned_list, bubble_in=.TRUE.)
  CALL ConstructMatrixFromTripletList(matAB, sorted_pruned_list, mat_c_rows, &
       & mat_c_columns)
  CALL DestructTripletList(sorted_pruned_list)
  CALL DestructTripletList(unsorted_pruned_list)
