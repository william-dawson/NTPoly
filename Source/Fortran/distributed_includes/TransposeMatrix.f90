  !! Local Data
  INTEGER :: II

  CALL ConstructTripletList(new_list)

  CALL GetMatrixTripletList(AMat, tlist)
  DO II = 1, tlist%CurrentSize
     IF (MOD(II, AMat%process_grid%num_process_slices) .EQ. &
          & AMat%process_grid%my_slice) THEN
        CALL GetTripletAt(tlist, II, trip)
        trip_t%index_row = trip%index_column
        trip_t%index_column = trip%index_row
        trip_t%point_value = trip%point_value
        CALL AppendToTripletList(new_list, trip_t)
     END IF
  END DO

  CALL DestructMatrix(TransMat)
  CALL ConstructEmptyMatrix(TransMat, AMat)
  CALL FillMatrixFromTripletList(TransMat, new_list)
  CALL DestructTripletList(new_list)
  CALL DestructTripletList(tlist)
