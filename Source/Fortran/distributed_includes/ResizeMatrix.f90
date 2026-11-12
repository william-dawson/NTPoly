  INTEGER :: II
  LOGICAL :: matrix_is_complex

  !! Get the triplet values.
  CALL GetMatrixTripletList(this, tlist)

  !! Prune the triplet values so that they fit in the new size.
  CALL ConstructTripletList(pruned)
  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, temp)
     IF (temp%index_row .LE. new_size .AND. &
          & temp%index_column .LE. new_size) THEN
        CALL AppendToTripletList(pruned, temp)
     END IF
  END DO

  !! Rebuild.
  matrix_is_complex = this%is_complex
  CALL ConstructEmptyMatrix(this, new_size, process_grid_in = this%process_grid, &
       & is_complex_in = matrix_is_complex)
  CALL FillMatrixFromTripletList(this, pruned, preduplicated_in = .TRUE.)

  !! Cleanup
  CALL DestructTripletList(tlist)
  CALL DestructTripletList(pruned)
