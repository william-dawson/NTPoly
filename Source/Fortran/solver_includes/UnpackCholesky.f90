  !! Local Variables
  INTEGER :: local_columns
  TYPE(TripletList_r) :: local_triplets
  TYPE(Triplet_r) :: temp
  INTEGER :: II, JJ

  local_columns = LMat%local_columns

  CALL ConstructTripletList(local_triplets)
  IF (LMat%process_grid%my_slice .EQ. 0) THEN
     DO JJ = 1, local_columns
        !! note transpose
        temp%index_row = JJ + LMat%start_column - 1
        DO II = 1, values_per_column(JJ)
           !! note transpose
           temp%index_column = INDEX(II, JJ) + LMat%start_row - 1
           temp%point_value = values(II, JJ)
           CALL AppendToTripletList(local_triplets, temp)
        END DO
     END DO
  END IF
  CALL FillMatrixFromTripletList(LMat, local_triplets)

  !! Cleanup
  CALL DestructTripletList(local_triplets)
