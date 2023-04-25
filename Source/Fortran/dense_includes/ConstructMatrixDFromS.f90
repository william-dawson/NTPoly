  !! Helper Variables
  INTEGER :: II, JJ
  INTEGER :: KK  ! Total element counter
  INTEGER :: elements_per_inner

  CALL ConstructEmptyMatrix(dense_matrix, sparse_matrix%rows, &
       & sparse_matrix%columns)
  dense_matrix%DATA = 0

  !! Loop over elements.
  KK = 1
  DO JJ = 1, sparse_matrix%columns
     elements_per_inner = sparse_matrix%outer_index(JJ + 1) - &
          & sparse_matrix%outer_index(JJ)
     temp%index_column = JJ
     DO II = 1, elements_per_inner
        temp%index_row = sparse_matrix%inner_index(KK)
        temp%point_value = sparse_matrix%values(KK)
        dense_matrix%DATA(temp%index_row, temp%index_column) = temp%point_value
        KK = KK + 1
     END DO
  END DO
