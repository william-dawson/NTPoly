  IF (.NOT. in%is_complex .EQV. convert_to_complex) THEN
     CALL MergeMatrixLocalBlocks(in, local_matrix)

     CALL ConstructEmptyMatrix(out, in%actual_matrix_dimension, &
          & process_grid_in=in%process_grid, is_complex_in=convert_to_complex)

     CALL ConvertMatrixType(local_matrix, converted_matrix)
     CALL SplitMatrixToLocalBlocks(out, converted_matrix)
  END IF
