  !! Local variables
  INTEGER :: number_of_values
  INTEGER :: start_index
  INTEGER :: II

  !! Allocate Memory
  CALL ConstructEmptyMatrix(column_out, this%rows, 1)
  start_index = this%outer_index(column_number)
  number_of_values = this%outer_index(column_number + 1) - &
       & this%outer_index(column_number)
  ALLOCATE(column_out%inner_index(number_of_values))
  ALLOCATE(column_out%values(number_of_values))

  !! Copy Values
  column_out%outer_index(1) = 0
  column_out%outer_index(2) = number_of_values
  DO II = 1, number_of_values
     column_out%inner_index(II) = this%inner_index(start_index + II)
     column_out%values(II) = this%values(start_index + II)
  END DO
