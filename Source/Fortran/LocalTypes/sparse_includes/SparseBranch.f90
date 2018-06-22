  !! Block A and B
  IF (.NOT. IsATransposed) THEN
     CALL matAT%Transpose(matA)
  END IF
  IF (.NOT. IsBTransposed) THEN
     CALL matBT%Transpose(matB)
  END IF

  IF (IsATransposed .AND. IsBTransposed) THEN
     CALL MultiplyBlock(matA, matB, blocked_memory_pool)
  ELSEIF (IsATransposed) THEN
     CALL MultiplyBlock(matA, matBT, blocked_memory_pool)
  ELSEIF (IsBTransposed) THEN
     CALL MultiplyBlock(matAT, matB, blocked_memory_pool)
  ELSE
     CALL MultiplyBlock(matAT, matBT, blocked_memory_pool)
  END IF

  !! Go from triplets to return matrix
  CALL PruneList(blocked_memory_pool, alpha, threshold, &
       & blocked_memory_pool%columns, blocked_memory_pool%rows, matC)
