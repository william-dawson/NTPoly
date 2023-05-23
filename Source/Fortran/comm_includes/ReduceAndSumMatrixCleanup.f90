  !! Local Data
  INTEGER :: II
  INTEGER :: total_values

  !! Build Matrix Objects
  CALL ConstructEmptyMatrix(acc_matrix, matrix%rows, matrix%columns)
  CALL ConstructEmptyMatrix(sum_matrix, matrix%rows, matrix%columns, &
       & zero_in = .TRUE.)

  !! Sum
  DO II = 1, helper%comm_size
     total_values = helper%values_per_process(II)
     ALLOCATE(acc_matrix%values(total_values))
     ALLOCATE(acc_matrix%inner_index(total_values))
     acc_matrix%values(:) = gathered_matrix%values( &
          & helper%displacement(II) + 1: &
          & helper%displacement(II) + helper%values_per_process(II))
     acc_matrix%inner_index(:) = gathered_matrix%inner_index( &
          & helper%displacement(II) + 1: &
          & helper%displacement(II) + helper%values_per_process(II))
     acc_matrix%outer_index(:) = gathered_matrix%outer_index(&
          & (matrix%columns + 1) * (II - 1) + 1:(matrix%columns + 1) * II)
     IF (II .EQ. helper%comm_size) THEN
        CALL IncrementMatrix(acc_matrix, sum_matrix, &
             & threshold_in = threshold)
     ELSE
        CALL IncrementMatrix(acc_matrix, sum_matrix,&
             & threshold_in = 0.0_NTREAL)
     END IF
     DEALLOCATE(acc_matrix%values)
     DEALLOCATE(acc_matrix%inner_index)
  END DO
  CALL CopyMatrix(sum_matrix, gathered_matrix)
  CALL DestructMatrix(sum_matrix)

  CALL DestructMatrix(acc_matrix)
  DEALLOCATE(helper%values_per_process)
  DEALLOCATE(helper%displacement)
