  !! Local Data
  INTEGER :: II, JJ
  INTEGER :: offset

  !! Sum Up The Outer Indices
  DO II = 1, helper%comm_size - 1
     offset = II * matrix%columns + 1
     DO JJ = 1, matrix%columns
        gathered_matrix%outer_index(offset + JJ) = &
             & gathered_matrix%outer_index(offset) + &
             & gathered_matrix%outer_index(offset + JJ)
     END DO
  END DO
  DEALLOCATE(helper%values_per_process)
  DEALLOCATE(helper%displacement)
