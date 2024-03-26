  REAL(NTREAL), PARAMETER :: sparsity_threshold = 0.1_NTREAL
  !! Counters and temporary data
  INTEGER :: mat_c_columns, mat_c_rows
  !! For Efficiency Purposes
  REAL(NTREAL) :: sparsity_a, sparsity_b
  REAL(NTREAL) :: sparsity_estimate
  LOGICAL :: pool_flag

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0_NTREAL
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(beta_in)) THEN
     beta = 0.0_NTREAL
  ELSE
     beta = beta_in
  END IF
  IF (.NOT. PRESENT(IsATransposed_in)) THEN
     IsATransposed = .FALSE.
  ELSE
     IsATransposed = IsATransposed_in
  END IF
  IF (.NOT. PRESENT(IsBTransposed_in)) THEN
     IsBTransposed = .FALSE.
  ELSE
     IsBTransposed = IsBTransposed_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0
  ELSE
     threshold = threshold_in
  END IF

  !! Storage details for result matrix
  IF (IsATransposed) THEN
     mat_c_rows = matA%columns
  ELSE
     mat_c_rows = matA%rows
  END IF
  IF (IsBTransposed) THEN
     mat_c_columns = matB%rows
  ELSE
     mat_c_columns = matB%columns
  END IF

  !! Initialization of Memory
  sparsity_a = DBLE(SIZE(matA%values)) / (matA%rows * matA%columns)
  sparsity_b = DBLE(SIZE(matB%values)) / (matB%rows * matB%columns)
  sparsity_estimate = 4*MAX(sparsity_a, sparsity_b)
  IF (sparsity_estimate > 1.0) THEN
     sparsity_estimate = 1.0
  ELSE IF (sparsity_estimate < 1e-8) THEN
     sparsity_estimate = 1e-8
  END IF

  CALL RoundTripLowP(matA, matAL)
  CALL RoundTripLowP(matB, matBL)

  !! Decide whether to do dense or sparse version.
  IF (MIN(sparsity_a, sparsity_b) .GT. sparsity_threshold) THEN
     CALL DenseBranch(matAL, matBL, matAB, IsATransposed, IsBTransposed, &
          & alpha, threshold)
  ELSE
     !! Setup the memory pool
     IF (.NOT. PRESENT(blocked_memory_pool_in)) THEN
        CALL ConstructMatrixMemoryPool(blocked_memory_pool, mat_c_columns, &
             & mat_c_rows, sparsity_estimate)
        pool_flag = .FALSE.
     ELSEIF (.NOT. CheckMemoryPoolValidity(blocked_memory_pool_in, &
          & mat_c_columns, mat_c_rows)) THEN
        CALL DestructMatrixMemoryPool(blocked_memory_pool_in)
        CALL ConstructMatrixMemoryPool(blocked_memory_pool_in, mat_c_columns, &
             & mat_c_rows, sparsity_estimate)
        pool_flag = .TRUE.
     ELSE
        CALL SetPoolSparsity(blocked_memory_pool_in, sparsity_estimate)
        pool_flag = .TRUE.
     END IF
     !! Multiply
     IF (pool_flag) THEN
        CALL SparseBranch(matAL, matBL, matAB, IsATransposed, IsBTransposed, &
             & alpha, threshold, blocked_memory_pool_in)
     ELSE
        CALL SparseBranch(matAL, matBL, matAB, IsATransposed, IsBTransposed, &
             & alpha, threshold, blocked_memory_pool)
     END IF
  END IF

  !! Handle the add part of GEMM
  IF (PRESENT(beta_in)) THEN
     IF (ABS(beta_in) .GT. 0) THEN
        CALL ScaleMatrix(matC, beta)
        CALL IncrementMatrix(matAB, matC)
     ELSE
        CALL CopyMatrix(matAB, matC)
     END IF
  ELSE
     CALL CopyMatrix(matAB, matC)
  END IF

  CALL DestructMatrix(matAB)
  CALL DestructMatrix(matAL)
  CALL DestructMatrix(matBL)
  CALL DestructMatrixMemoryPool(blocked_memory_pool)
