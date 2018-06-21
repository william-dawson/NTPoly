  !! Allocate New Matrix
  num_values = matA%outer_index(matA%columns+1)
  CALL this%InitEmpty(matA%columns, matA%rows)
  ALLOCATE(this%inner_index(num_values))
  ALLOCATE(this%values(num_values))

  !! Temporary Arrays
  ALLOCATE(values_per_row(matA%rows))
  ALLOCATE(offset_array(matA%rows))

  !! Count the values per row
  values_per_row = 0
  DO II = 1, num_values
     inner_index = matA%inner_index(II)
     values_per_row(inner_index) = values_per_row(inner_index) + 1
  END DO
  offset_array(1) = 0
  DO II = 2, matA%rows
     offset_array(II) = offset_array(II-1) + values_per_row(II-1)
  END DO

  !! Insert
  this%outer_index(:matA%rows) = offset_array(:matA%rows)
  this%outer_index(matA%rows+1) = offset_array(matA%rows) + &
       & values_per_row(matA%rows)
  DO II = 1, matA%columns
     elements_per_inner = matA%outer_index(II+1) - matA%outer_index(II)
     matA_offset = matA%outer_index(II)
     DO JJ = 1, elements_per_inner
        inner_index = matA%inner_index(matA_offset+JJ)
        insert_pt = offset_array(inner_index)+1
        this%inner_index(insert_pt) = II
        this%values(insert_pt) = matA%values(matA_offset+JJ)
        offset_array(inner_index) = offset_array(inner_index) +1
     END DO
  END DO

  !! Cleanup
  DEALLOCATE(values_per_row)
  DEALLOCATE(offset_array)
