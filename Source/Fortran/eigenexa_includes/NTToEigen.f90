  !! Local Variables
  INTEGER :: ilookup, jlookup, iowner, jowner, ijowner
  INTEGER :: lrow, lcol
  INTEGER :: II

  !! We will fill a triplet list for each other process
  ALLOCATE(send_list(exa%num_procs))
  DO II = 1, exa%num_procs
     CALL ConstructTripletList(send_list(II))
  END DO

  !! Now Get The Triplet List, and adjust
  CALL GetMatrixTripletList(A, triplet_a)
  DO II = 1, triplet_a%CurrentSize
     CALL GetTripletAt(triplet_a, II, trip)

     !! Determine where that triplet will reside
     iowner = eigen_owner_node(trip%index_row, exa%proc_rows, exa%rowid)
     jowner = eigen_owner_node(trip%index_column, exa%proc_cols, exa%colid)
     ijowner = (jowner-1)*exa%proc_rows + iowner

     !! New indices
     ilookup = eigen_translate_g2l(trip%index_row, exa%proc_rows, &
          & exa%rowid)
     jlookup = eigen_translate_g2l(trip%index_column, exa%proc_cols, &
          & exa%colid)
     CALL SetTriplet(shifted_trip, jlookup, ilookup, trip%point_value)

     CALL AppendToTripletList(send_list(ijowner), shifted_trip)
  END DO

  !! Redistribute The Triplets
  CALL RedistributeTripletLists(send_list, exa%comm, recv_list)

  !! Write To The Dense Array
  DO II = 1, recv_list%CurrentSize
     lrow = recv_list%DATA(II)%index_row
     lcol = recv_list%DATA(II)%index_column
     AD(lrow,lcol) = recv_list%DATA(II)%point_value
  END DO

  !! Cleanup
  DO II = 1, exa%num_procs
     CALL DestructTripletList(send_list(II))
  END DO
  DEALLOCATE(send_list)
  CALL DestructTripletList(recv_list)
  CALL DestructTripletList(triplet_a)
