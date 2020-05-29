  !! Counter Variables
  INTEGER :: outer_counter
  INTEGER :: inner_a, inner_b
  INTEGER :: total_counter_a, total_counter_b, total_counter_c
  !! Temporary Variables
  INTEGER :: indices_added_into_c
  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  INTEGER :: size_of_a, size_of_b

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0_NTREAL
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0_NTREAL
  ELSE
     threshold = threshold_in
  END IF

  size_of_a = matA%outer_index(matA%columns+1)

  !! Allocate sufficient space for matC
  CALL ConstructEmptyMatrix(matC, matA%rows, matA%columns)
  IF (ALLOCATED(matB%values)) THEN
     size_of_b = matB%outer_index(matB%columns+1)
     ALLOCATE(matC%inner_index(size_of_a+size_of_b))
     ALLOCATE(matC%values(size_of_a+size_of_b))
  ELSE
     ALLOCATE(matC%inner_index(size_of_a))
     ALLOCATE(matC%values(size_of_a))
  END IF

  !! Perform loops
  total_counter_a = 1
  total_counter_b = 1
  total_counter_c = 1
  DO outer_counter = 1, matA%columns
     !! Inner counters
     inner_a = matA%outer_index(outer_counter+1) - &
          & matA%outer_index(outer_counter)
     inner_b = matB%outer_index(outer_counter+1) - &
          & matB%outer_index(outer_counter)
     CALL AddSparseVectors(&
          matA%inner_index(total_counter_a:total_counter_a+inner_a-1),&
          matA%values(total_counter_a:total_counter_a+inner_a-1),&
          matB%inner_index(total_counter_b:total_counter_b+inner_b-1),&
          matB%values(total_counter_b:total_counter_b+inner_b-1),&
          matC%inner_index(total_counter_c:),matC%values(total_counter_c:),&
          indices_added_into_c, alpha, threshold)
     matC%outer_index(outer_counter+1) = matC%outer_index(outer_counter)+&
          & indices_added_into_c
     total_counter_a = total_counter_a + inner_a
     total_counter_b = total_counter_b + inner_b
     total_counter_c = total_counter_c + indices_added_into_c
  END DO

  !! Cleanup
  CALL DestructMatrix(matB)
  CALL ConstructEmptyMatrix(matB, matC%rows, matC%columns)
  matB%outer_index = matC%outer_index
  ALLOCATE(matB%inner_index(matC%outer_index(matC%columns+1)))
  ALLOCATE(matB%values(matC%outer_index(matC%columns+1)))
  matB%inner_index = matC%inner_index(:matC%outer_index(matC%columns+1))
  matB%values = matC%values(:matC%outer_index(matC%columns+1))
  CALL DestructMatrix(matC)
