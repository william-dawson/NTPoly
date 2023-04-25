  !! Local Data
  INTEGER :: II
  INTEGER :: temporary_total_values

  !! Build Matrix Objects
  CALL ConstructEmptyMatrix(temporary_matrix, matrix%rows, matrix%columns)
  CALL ConstructEmptyMatrix(sum_matrix, matrix%rows, matrix%columns, &
       & zero_in=.TRUE.)

  !! Sum
  DO II = 1, helper%comm_size
     temporary_total_values = helper%values_per_process(II)
     ALLOCATE(temporary_matrix%values(temporary_total_values))
     ALLOCATE(temporary_matrix%inner_index(temporary_total_values))
     temporary_matrix%values = gathered_matrix%values( &
          & helper%displacement(II) + 1: &
          & helper%displacement(II) + helper%values_per_process(II))
     temporary_matrix%inner_index = gathered_matrix%inner_index( &
          & helper%displacement(II) + 1: &
          & helper%displacement(II) + helper%values_per_process(II))
     temporary_matrix%outer_index = gathered_matrix%outer_index(&
          & (matrix%columns + 1)*(II - 1) + 1:(matrix%columns + 1)*(II))
     IF (II .EQ. helper%comm_size) THEN
        CALL IncrementMatrix(temporary_matrix, sum_matrix, &
             & threshold_in=threshold)
     ELSE
        CALL IncrementMatrix(temporary_matrix,sum_matrix,&
             & threshold_in=REAL(0.0, NTREAL))
     END IF
     DEALLOCATE(temporary_matrix%values)
     DEALLOCATE(temporary_matrix%inner_index)
  END DO
  CALL CopyMatrix(sum_matrix, gathered_matrix)
  CALL DestructMatrix(sum_matrix)

  CALL DestructMatrix(temporary_matrix)
  DEALLOCATE(helper%values_per_process)
  DEALLOCATE(helper%displacement)
