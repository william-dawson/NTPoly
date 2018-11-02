
  INTEGER :: II

  CALL ConstructTripletList(outlist)
  DO II = 1, inlist%CurrentSize
     CALL GetTripletAt(inlist, II, temp)
     CALL proc(temp%index_row, temp%index_column, temp%point_value)
     CALL AppendToTripletList(outlist, temp)
  END DO
