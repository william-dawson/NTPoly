  TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
  !! Local Data
  INTEGER :: counter, inner_counter
  INTEGER :: temp_offset

  !! Sum Up The Outer Indices
  DO counter = 1, helper%comm_size - 1
     temp_offset = counter*matrix%columns+1
     DO inner_counter = 1, matrix%columns
        gathered_matrix%outer_index(temp_offset+inner_counter) = &
             & gathered_matrix%outer_index(temp_offset) + &
             & gathered_matrix%outer_index(temp_offset+inner_counter)
     END DO
  END DO
  DEALLOCATE(helper%values_per_process)
  DEALLOCATE(helper%displacement)
