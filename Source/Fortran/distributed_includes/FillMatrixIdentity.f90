  !! Local Data
  INTEGER :: i, j
  INTEGER :: total_values

  !! There can't be more than one entry per row
  CALL ConstructTripletList(triplet_list, this%local_rows)

  total_values = 0
  !! Find local identity values
  row_iter: DO j = 1, this%local_rows
     column_iter: DO i = 1, this%local_columns
        IF (j + this%start_row - 1 .EQ. i + this%start_column - 1 .AND. &
             & j+this%start_row-1 .LE. this%actual_matrix_dimension) THEN
           total_values = total_values + 1
           triplet_list%data(total_values)%index_column = i + this%start_column - 1 
           triplet_list%data(total_values)%index_row = j + this%start_row - 1
           triplet_list%data(total_values)%point_value = 1.0
        END IF
     END DO column_iter
  END DO row_iter
  triplet_list%CurrentSize = total_values

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, unsorted_triplet_list, prepartitioned_in=.TRUE.)

  CALL DestructTripletList(triplet_list)
