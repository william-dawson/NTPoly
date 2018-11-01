  INTEGER :: II

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
  CALL ConstructEmptyMatrix(this, new_size)
  CALL FillMatrixFromTripletList(this, pruned, preduplicated_in=.TRUE.)

  !! Cleanup
  CALL DestructTripletList(tlist)
  CALL DestructTripletList(pruned)
