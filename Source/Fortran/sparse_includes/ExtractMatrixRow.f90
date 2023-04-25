  !! Temporary Variables
  INTEGER :: values_found
  INTEGER :: elements_per_inner
  INTEGER :: II, JJ, KK

  !! Fill a value buffer
  CALL ConstructEmptyMatrix(row_out, 1, this%columns)
  ALLOCATE(value_buffer(this%columns))
  values_found = 0
  KK = 1
  row_out%outer_index(1) = 0
  DO II = 1, this%columns
     row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + &
          & row_out%outer_index(II)
     elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
     DO JJ = 1, elements_per_inner
        IF (this%inner_index(KK) .EQ. row_number) THEN
           values_found = values_found + 1
           value_buffer(values_found) = this%values(KK)
           row_out%outer_index(II + 1) = row_out%outer_index(II + 1) + 1
        END IF
        KK = KK + 1
     END DO
  END DO

  !! Copy To Actual Matrix
  ALLOCATE(row_out%inner_index(values_found))
  row_out%inner_index = 1
  ALLOCATE(row_out%values(values_found))
  row_out%values = value_buffer(:values_found)

  !! Cleanup
  DEALLOCATE(value_buffer)
