  !! Local Data
  INTEGER :: counter
  INTEGER :: size_temp
  TYPE(ProcessGrid_t) :: grid_temp
  LOGICAL :: is_complex_temp

  CALL GetMatrixTripletList(this, triplet_list)
  CALL ConstructTripletList(new_list)

  DO counter=1,triplet_list%CurrentSize
     CALL GetTripletAt(triplet_list, counter, temporary)
     IF (ABS(temporary%point_value) .GT. threshold) THEN
        CALL AppendToTripletList(new_list, temporary)
     END IF
  END DO

  size_temp = this%actual_matrix_dimension
  grid_temp = this%process_grid
  is_complex_temp = this%is_complex
  CALL DestructMatrix(this)
  CALL ConstructEmptyMatrix(this, size_temp, grid_temp, is_complex_temp)

  CALL FillMatrixFromTripletList(this, new_list)
