  !! Local Data
  INTEGER :: counter
  INTEGER :: temporary_total_values

  !! Build Matrix Objects
  CALL temporary_matrix%InitEmpty(matrix%rows,matrix%columns)
  CALL temporary_matrix%InitEmpty(matrix%rows,matrix%columns)

  !! Sum
  DO counter = 1, helper%comm_size
     temporary_total_values = helper%values_per_process(counter)
     ALLOCATE(temporary_matrix%values(temporary_total_values))
     ALLOCATE(temporary_matrix%inner_index(temporary_total_values))
     temporary_matrix%values = gathered_matrix%values( &
          & helper%displacement(counter)+1: &
          & helper%displacement(counter) + helper%values_per_process(counter))
     temporary_matrix%inner_index = gathered_matrix%inner_index( &
          & helper%displacement(counter)+1: &
          & helper%displacement(counter) + helper%values_per_process(counter))
     temporary_matrix%outer_index = gathered_matrix%outer_index(&
          & (matrix%columns+1)*(counter-1)+1:(matrix%columns+1)*(counter))
     IF (counter .EQ. helper%comm_size) THEN
        CALL IncrementMatrix(temporary_matrix, sum_matrix, &
             & threshold_in=threshold)
     ELSE
        CALL IncrementMatrix(temporary_matrix, sum_matrix)
     END IF
     DEALLOCATE(temporary_matrix%values)
     DEALLOCATE(temporary_matrix%inner_index)
  END DO
  CALL gathered_matrix%Copy(sum_matrix)
  CALL sum_matrix%Destruct

  CALL temporary_Matrix%Destruct
  DEALLOCATE(helper%values_per_process)
  DEALLOCATE(helper%displacement)
