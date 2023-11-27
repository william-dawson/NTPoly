  !! Local Data
  INTEGER :: II, KK

  CALL GetMatrixTripletList(AMat, tlist)

  !! Determine the size.
  !! these loops could be merged using append but I choose not to for the
  !! sake of memory.
  KK = 0
  DO II = 1, tlist%CurrentSize
     IF (MOD(II, AMat%process_grid%num_process_slices) .EQ. &
          & AMat%process_grid%my_slice) THEN
        KK = KK + 1
     END IF
  END DO
  CALL ConstructTripletList(new_list, KK)

  !! Fill the Triplet List
  KK = 0
  DO II = 1, tlist%CurrentSize
     IF (MOD(II, AMat%process_grid%num_process_slices) .EQ. &
          & AMat%process_grid%my_slice) THEN
        KK = KK + 1
        CALL GetTripletAt(tlist, II, trip)
        new_list%data(KK)%index_row = trip%index_column
        new_list%data(KK)%index_column = trip%index_row
        new_list%data(KK)%point_value = trip%point_value
     END IF
  END DO
  CALL DestructTripletList(tlist)

  !! Create the matrix
  CALL DestructMatrix(TransMat)
  CALL ConstructEmptyMatrix(TransMat, AMat)
  CALL FillMatrixFromTripletList(TransMat, new_list)
  CALL DestructTripletList(new_list)