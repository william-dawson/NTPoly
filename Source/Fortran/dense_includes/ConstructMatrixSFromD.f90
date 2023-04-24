  !! Local Variables
  INTEGER :: II, JJ, KK, NNZ
  REAL(NTREAL) :: threshold

  IF (PRESENT(threshold_in)) THEN
     threshold = threshold_in
  ELSE
     threshold = 0.0_NTREAL
  END IF

  CALL ConstructEmptyMatrix(sparse_matrix, dense_matrix%rows, &
       & dense_matrix%columns)

  !! Fill in the outer index information.
  NNZ = 0
  DO II = 1, dense_matrix%columns
     DO JJ = 1, dense_matrix%rows
        IF (ABS(dense_matrix%DATA(JJ, II)) .GT. threshold) THEN
           NNZ = NNZ + 1
        END IF
     END DO
     sparse_matrix%outer_index(II + 1) = NNZ
  END DO

  !! Allocate Storage
  ALLOCATE(sparse_matrix%inner_index(NNZ))
  ALLOCATE(sparse_matrix%values(NNZ))

  !! Fill in the Values
  KK = 1
  DO II = 1, dense_matrix%columns
     DO JJ = 1, dense_matrix%rows
        IF (ABS(dense_matrix%DATA(JJ, II)) .GT. threshold) THEN
           sparse_matrix%inner_index(KK) = JJ
           sparse_matrix%values(KK) = dense_matrix%DATA(JJ, II)
           KK = KK + 1
        END IF
     END DO
  END DO
