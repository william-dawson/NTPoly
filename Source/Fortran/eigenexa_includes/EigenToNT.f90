  !! Local Variables
  INTEGER :: row_start, row_end, col_start, col_end
  INTEGER :: II, JJ, ilookup, jlookup
  INTEGER :: ind

  !! Get The Eigenvectors
  row_start = eigen_loop_start(1, exa%proc_rows, exa%rowid)
  row_end = eigen_loop_end(exa%mat_dim, exa%proc_rows, exa%rowid)
  col_start = eigen_loop_start(1, exa%proc_cols, exa%colid)
  col_end = eigen_loop_end(exa%mat_dim, exa%proc_cols, exa%colid)

  !! Convert to a 1D array for index ease.
  ALLOCATE(VD1(SIZE(VD, DIM = 1)*SIZE(VD, DIM = 2)))
  VD1 = PACK(VD, .TRUE.)

  CALL ConstructTripletList(triplet_v)
  ind = 1
  DO JJ = col_start, col_end
     jlookup = eigen_translate_l2g(JJ, exa%proc_cols, exa%colid)
     DO II = row_start, row_end
        IF (ABS(VD1(ind + II -1)) .GT. params%threshold) THEN
           ilookup = eigen_translate_l2g(II, exa%proc_rows, exa%rowid)
           CALL SetTriplet(trip, jlookup, ilookup, VD1(ind + II -1))
           CALL AppendToTripletList(triplet_v, trip)
        END IF
     END DO
     ind = ind + exa%offset
  END DO

  CALL FillMatrixFromTripletList(V, triplet_v)

  !! Cleanup
  CALL DestructTripletList(triplet_v)

  DEALLOCATE(VD1)
