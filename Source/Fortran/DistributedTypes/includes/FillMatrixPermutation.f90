  !! Local Data
  INTEGER :: total_values
  INTEGER :: counter
  INTEGER :: local_row, local_column

  !! Build Local Triplet List
  !! There can't be more than one entry per row
  CALL ConstructTripletList(triplet_list, this%local_rows)
  total_values = 0
  IF (rows) THEN
     DO counter=this%start_row,this%end_row-1
        IF (permutation_vector(counter) .GE. this%start_column .AND. &
             & permutation_vector(counter) .LT. this%end_column) THEN
           total_values = total_values + 1
           local_column = permutation_vector(counter) - this%start_column + 1
           local_row = counter - this%start_row + 1
           triplet_list%data(total_values)%index_column = local_column
           triplet_list%data(total_values)%index_row = local_row
           triplet_list%data(total_values)%point_value = 1.0
        END IF
     END DO
  ELSE
     DO counter=this%start_column,this%end_column-1
        IF (permutation_vector(counter) .GE. this%start_row .AND. &
             & permutation_vector(counter) .LT. this%end_row) THEN
           total_values = total_values + 1
           local_column = counter - this%start_column + 1
           local_row = permutation_vector(counter) - this%start_row + 1
           triplet_list%data(total_values)%index_column = local_column
           triplet_list%data(total_values)%index_row = local_row
           triplet_list%data(total_values)%point_value = 1.0
        END IF
     END DO
  END IF

  !! Finish constructing
  CALL ConstructTripletList(unsorted_triplet_list, total_values)
  unsorted_triplet_list%data = triplet_list%data(:total_values)
  CALL SortTripletList(unsorted_triplet_list, this%local_columns, &
       & this%local_rows, sorted_triplet_list)
  CALL ConstructMatrixFromTripletList(local_matrix, sorted_triplet_list, &
       & this%local_rows, this%local_columns)

  CALL SplitMatrixToLocalBlocks(this, local_matrix)

  CALL DestructMatrix(local_matrix)
  CALL DestructTripletList(triplet_list)
  CALL DestructTripletList(unsorted_triplet_list)
  CALL DestructTripletList(sorted_triplet_list)
