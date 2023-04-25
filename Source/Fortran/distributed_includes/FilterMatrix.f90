  !! Local Data
  TYPE(ProcessGrid_t) :: grid_temp
  LOGICAL :: is_complex_temp
  INTEGER :: II, size_temp

  CALL GetMatrixTripletList(this, tlist)
  CALL ConstructTripletList(new_list)

  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, trip)
     IF (ABS(trip%point_value) .GT. threshold) THEN
        CALL AppendToTripletList(new_list, trip)
     END IF
  END DO

  size_temp = this%actual_matrix_dimension
  grid_temp = this%process_grid
  is_complex_temp = this%is_complex
  CALL DestructMatrix(this)
  CALL ConstructEmptyMatrix(this, size_temp, grid_temp, is_complex_temp)

  CALL FillMatrixFromTripletList(this, new_list)
