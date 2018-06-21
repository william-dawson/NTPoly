  !! Allocate New Matrix
  num_values = this%outer_index(this%columns+1)
  CALL matT%InitEmpty(this%columns, this%rows)
  ALLOCATE(matT%inner_index(num_values))
  ALLOCATE(matT%values(num_values))

  !! Temporary Arrays
  ALLOCATE(values_per_row(this%rows))
  ALLOCATE(offset_array(this%rows))

  !! Count the values per row
  values_per_row = 0
  DO II = 1, num_values
     inner_index = this%inner_index(II)
     values_per_row(inner_index) = values_per_row(inner_index) + 1
  END DO
  offset_array(1) = 0
  DO II = 2, this%rows
     offset_array(II) = offset_array(II-1) + values_per_row(II-1)
  END DO

  !! Insert
  matT%outer_index(:this%rows) = offset_array(:this%rows)
  matT%outer_index(this%rows+1) = offset_array(this%rows) + &
       & values_per_row(this%rows)
  DO II = 1, this%columns
     elements_per_inner = this%outer_index(II+1) - this%outer_index(II)
     this_offset = this%outer_index(II)
     DO JJ = 1, elements_per_inner
        inner_index = this%inner_index(this_offset+JJ)
        insert_pt = offset_array(inner_index)+1
        matT%inner_index(insert_pt) = II
        matT%values(insert_pt) = this%values(this_offset+JJ)
        offset_array(inner_index) = offset_array(inner_index) +1
     END DO
  END DO

  !! Cleanup
  DEALLOCATE(values_per_row)
  DEALLOCATE(offset_array)
