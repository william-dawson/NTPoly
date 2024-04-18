  !! Local Data
  INTEGER :: II

  CALL GetMatrixTripletList(this, tlist)
  CALL ConstructTripletList(new_list)

  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, trip)
     IF (ABS(trip%point_value) .GT. threshold) THEN
        CALL AppendToTripletList(new_list, trip)
     END IF
  END DO

  CALL FillMatrixFromTripletList(this, new_list, preduplicated_in=.TRUE.)
