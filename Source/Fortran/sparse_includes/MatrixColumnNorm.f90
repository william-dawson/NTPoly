  !! Local Data
  INTEGER :: outer_counter, inner_counter
  INTEGER :: elements_per_inner

  !! Allocate Space For Result
  norm_per_column = 0

  !! Iterate Over Local Data
  DO outer_counter = 1, this%columns
     elements_per_inner = this%outer_index(outer_counter + 1) - &
          & this%outer_index(outer_counter)
     DO inner_counter = 1, elements_per_inner
        temp_value = this%values(this%outer_index(outer_counter)+ &
             & inner_counter)
        norm_per_column(outer_counter) = norm_per_column(outer_counter) + &
             & ABS(temp_value)
     END DO
  END DO
