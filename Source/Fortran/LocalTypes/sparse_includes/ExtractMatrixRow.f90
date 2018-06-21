  !! Fill a value buffer
  CALL row_out%InitEmpty(1, this%columns)
  ALLOCATE(value_buffer(this%columns))
  values_found = 0
  total_counter = 1
  row_out%outer_index(1) = 0
  DO outer_counter = 1, this%columns
     row_out%outer_index(outer_counter+1) = &
          & row_out%outer_index(outer_counter+1) + &
          & row_out%outer_index(outer_counter)
     elements_per_inner = this%outer_index(outer_counter+1) - &
          & this%outer_index(outer_counter)
     DO inner_counter = 1, elements_per_inner
        IF (this%inner_index(total_counter) .EQ. row_number) THEN
           values_found = values_found + 1
           value_buffer(values_found) = this%values(total_counter)
           row_out%outer_index(outer_counter+1) = &
                & row_out%outer_index(outer_counter+1) + 1
        END IF
        total_counter = total_counter + 1
     END DO
  END DO

  !! Copy To Actual Matrix
  ALLOCATE(row_out%inner_index(values_found))
  row_out%inner_index = 1
  ALLOCATE(row_out%values(values_found))
  row_out%values = value_buffer(:values_found)

  !! Cleanup
  DEALLOCATE(value_buffer)
