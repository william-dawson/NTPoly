  !! Counter Variables
  INTEGER :: II
  INTEGER :: inner_a, inner_b
  INTEGER :: total_a, total_b, total_c
  !! Temporary Variables
  INTEGER :: indices_added_into_c
  INTEGER :: size_of_a, size_of_b

  CALL ConstructEmptyMatrix(TempMat, matA%rows, matA%columns)
  size_of_a = matA%outer_index(matA%columns + 1)
  size_of_b = matB%outer_index(matB%columns + 1)
  ALLOCATE(TempMat%inner_index(MIN(size_of_a, size_of_b)))
  ALLOCATE(TempMat%values(MIN(size_of_a, size_of_b)))

  !! Perform loops
  total_a = 1
  total_b = 1
  total_c = 1
  DO II = 1, matA%columns
     !! Inner counters
     inner_a = matA%outer_index(II + 1) - matA%outer_index(II)
     inner_b = matB%outer_index(II + 1) - matB%outer_index(II)
     CALL PairwiseMultiplyVectors(&
          matA%inner_index(total_a:total_a + inner_a - 1), &
          matA%values(total_a:total_a + inner_a - 1), &
          matB%inner_index(total_b:total_b + inner_b - 1), &
          matB%values(total_b:total_b + inner_b - 1), &
          TempMat%inner_index(total_c:), &
          TempMat%values(total_c:), &
          indices_added_into_c)
     TempMat%outer_index(II + 1) = &
          & TempMat%outer_index(II) + indices_added_into_c
     total_a = total_a + inner_a
     total_b = total_b + inner_b
     total_c = total_c + indices_added_into_c
  END DO

  !! Cleanup
  CALL DestructMatrix(matC)
  CALL ConstructEmptyMatrix(matC, TempMat%rows, TempMat%columns)
  matC%outer_index(:) = TempMat%outer_index
  ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
  ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
  matC%inner_index(:) = TempMat%inner_index(&
       & :TempMat%outer_index(TempMat%columns+1))
  matC%values(:) = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
  CALL DestructMatrix(TempMat)
