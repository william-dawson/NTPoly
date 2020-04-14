  !! Local Data
  TYPE(ReduceHelper_t) :: row_helper
  TYPE(ReduceHelper_t) :: column_helper

  CALL MergeMatrixLocalBlocks(this, local)

  !! Merge Columns
  CALL TransposeMatrix(local, localT)
  CALL ReduceAndComposeMatrixSizes(localT, this%process_grid%column_comm, &
       & merged_columns, column_helper)
  DO WHILE(.NOT. TestReduceSizeRequest(column_helper))
  END DO
  CALL ReduceAndComposeMatrixData(localT, this%process_grid%column_comm, &
       & merged_columns, column_helper)
  DO WHILE(.NOT. TestReduceInnerRequest(column_helper))
  END DO
  DO WHILE(.NOT. TestReduceDataRequest(column_helper))
  END DO
  CALL ReduceAndComposeMatrixCleanup(localT, merged_columns, column_helper)

  !! Merge Rows
  CALL TransposeMatrix(merged_columns, merged_columnsT)
  CALL ReduceAndComposeMatrixSizes(merged_columnsT, &
       & this%process_grid%row_comm, gathered, row_helper)
  DO WHILE(.NOT. TestReduceSizeRequest(row_helper))
  END DO
  CALL ReduceAndComposeMatrixData(merged_columnsT, &
       & this%process_grid%row_comm, gathered, row_helper)
  DO WHILE(.NOT. TestReduceInnerRequest(row_helper))
  END DO
  DO WHILE(.NOT. TestReduceDataRequest(row_helper))
  END DO
  CALL ReduceAndComposeMatrixCleanup(merged_columnsT, gathered, &
        & row_helper)

  !! Remove the excess rows and columns that come from the logical size.
  CALL ConstructEmptyMatrix(local_mat, this%actual_matrix_dimension, &
       & this%actual_matrix_dimension)
  local_mat%outer_index = gathered%outer_index(:this%actual_matrix_dimension+1)
  ALLOCATE(local_mat%inner_index(SIZE(gathered%inner_index)))
  local_mat%inner_index = gathered%inner_index
  ALLOCATE(local_mat%values(SIZE(gathered%values)))
  local_mat%values = gathered%values
