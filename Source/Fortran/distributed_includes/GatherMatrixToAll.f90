  CALL MergeMatrixLocalBlocks(this, local)

  !! Merge Columns
  CALL TransposeMatrix(local, localT)
  CALL ReduceAndComposeMatrix(localT, this%process_grid%column_comm, &
       & merged_columns)

  !! Merge Rows
  CALL TransposeMatrix(merged_columns, merged_columnsT)
  CALL ReduceAndComposeMatrix(merged_columnsT, this%process_grid%row_comm, &
       & gathered)

  !! Remove the excess rows and columns that come from the logical size.
  CALL ConstructEmptyMatrix(local_mat, this%actual_matrix_dimension, &
       & this%actual_matrix_dimension)
  local_mat%outer_index(:) = &
       & gathered%outer_index(:this%actual_matrix_dimension + 1)
  ALLOCATE(local_mat%inner_index(SIZE(gathered%inner_index)))
  local_mat%inner_index(:) = gathered%inner_index
  ALLOCATE(local_mat%values(SIZE(gathered%values)))
  local_mat%values(:) = gathered%values
