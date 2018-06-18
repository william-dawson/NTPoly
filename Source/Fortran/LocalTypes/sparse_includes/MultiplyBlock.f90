  INTEGER :: temp_inserted_values
  INTEGER :: temp_index_a, temp_index_b
  INTEGER :: elements_per_inner_a
  INTEGER :: elements_per_inner_b
  LOGICAL :: is_set
  !! Counters
  INTEGER :: outer_counter, inner_counter_a, inner_counter_b

  !! Multiply
  DO outer_counter = 1, matAT%columns
     elements_per_inner_a = matAT%outer_index(outer_counter+1) - &
          & matAT%outer_index(outer_counter)
     DO inner_counter_a = 1, elements_per_inner_a
        temp_value_a = matAT%values(matAT%outer_index(outer_counter)+ &
             & inner_counter_a)
        temp_index_a = matAT%inner_index(matAT%outer_index(outer_counter)+ &
             & inner_counter_a)
        elements_per_inner_b = matBT%outer_index(temp_index_a+1) - &
             & matBT%outer_index(temp_index_a)
        DO inner_counter_b = 1, elements_per_inner_b
           temp_index_b = matBT%inner_index(matBT%outer_index(temp_index_a)+ &
                & inner_counter_b)
           temp_value_b = matBT%values(matBT%outer_index(temp_index_a)+ &
                & inner_counter_b)
           temp_value_c = memorypool%value_array(temp_index_b,outer_counter)
           is_set = memorypool%dirty_array(temp_index_b,outer_counter)
           IF (is_set .EQV. .FALSE.) THEN
              memorypool%dirty_array(temp_index_b,outer_counter) = .TRUE.
              temp_inserted_values = memorypool%inserted_per_bucket(&
                   & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) + 1
              memorypool%inserted_per_bucket(&
                   & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) = &
                   & temp_inserted_values
              memorypool%hash_index(temp_inserted_values+ &
                   & ((temp_index_b-1)/memorypool%hash_size)*memorypool%hash_size, &
                   & outer_counter) = temp_index_b
           END IF
           memorypool%value_array(temp_index_b,outer_counter) = &
                & temp_value_c + temp_value_a*temp_value_b
        END DO
     END DO
  END DO
