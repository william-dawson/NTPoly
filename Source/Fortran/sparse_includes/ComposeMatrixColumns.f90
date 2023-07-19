  !! Local Variables
  INTEGER :: total_columns, total_values
  INTEGER :: inner_start, inner_length
  INTEGER :: outer_start, outer_length
  INTEGER :: outer_offset
  INTEGER :: II, JJ
  INTEGER :: size_of_mat

  CALL DestructMatrix(out_matrix)

  !! Figure Out The Sizes
  total_columns = 0
  total_values  = 0
  DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
     total_columns = total_columns + mat_list(II)%columns
     size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
     total_values  = total_values + size_of_mat
  END DO

  !! Allocate The Space
  CALL ConstructEmptyMatrix(out_matrix, mat_list(1)%rows, total_columns)
  ALLOCATE(out_matrix%inner_index(total_values))
  ALLOCATE(out_matrix%values(total_values))

  !! Fill In The Values
  inner_start = 1
  outer_start = 1
  outer_offset = 0
  DO II = LBOUND(mat_list, dim = 1), UBOUND(mat_list, dim = 1)
     !! Inner indices and values
     size_of_mat = mat_list(II)%outer_index(mat_list(II)%columns + 1)
     inner_length = size_of_mat
     DO JJ = 1, SIZE(mat_list(II)%inner_index)
       out_matrix%inner_index(inner_start + JJ - 1) = mat_list(II)%inner_index(JJ)
       out_matrix%values(inner_start + JJ - 1) = mat_list(II)%values(JJ)
     END DO
     inner_start = inner_start + inner_length
     !! Outer Indices
     outer_length = mat_list(II)%columns + 1
     out_matrix%outer_index(outer_start:outer_start+outer_length - 1) = &
          & mat_list(II)%outer_index + outer_offset
     outer_start = outer_start + outer_length - 1
     outer_offset = out_matrix%outer_index(outer_start)
  END DO
