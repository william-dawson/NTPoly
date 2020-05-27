  !! Counter Variables
  INTEGER :: outer_counter
  INTEGER :: inner_a, inner_b
  INTEGER :: total_counter_a, total_counter_b, total_counter_c
  !! Temporary Variables
  INTEGER :: indices_added_into_c
  INTEGER :: size_of_a, size_of_b

  CALL ConstructEmptyMatrix(TempMat, matA%rows, matA%columns)
  size_of_a = matA%outer_index(matA%columns+1)
  size_of_b = matB%outer_index(matB%columns+1)
  ALLOCATE(TempMat%inner_index(MIN(size_of_a,size_of_b)))
  ALLOCATE(TempMat%values(MIN(size_of_a,size_of_b)))

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
     CALL PairwiseMultiplyVectors(&
          matA%inner_index(total_counter_a:total_counter_a+inner_a-1),&
          matA%values(total_counter_a:total_counter_a+inner_a-1),&
          matB%inner_index(total_counter_b:total_counter_b+inner_b-1),&
          matB%values(total_counter_b:total_counter_b+inner_b-1),&
          TempMat%inner_index(total_counter_c:),&
          TempMat%values(total_counter_c:),&
          indices_added_into_c)
     TempMat%outer_index(outer_counter+1) = TempMat%outer_index(outer_counter)+&
          & indices_added_into_c
     total_counter_a = total_counter_a + inner_a
     total_counter_b = total_counter_b + inner_b
     total_counter_c = total_counter_c + indices_added_into_c
  END DO

  !! Cleanup
  CALL DestructMatrix(matC)
  CALL ConstructEmptyMatrix(matC, TempMat%rows, TempMat%columns)
  matC%outer_index = TempMat%outer_index
  ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
  ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
  matC%inner_index = TempMat%inner_index(&
       & :TempMat%outer_index(TempMat%columns+1))
  matC%values = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
  CALL DestructMatrix(TempMat)
