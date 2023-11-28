  !! Local Variables
  INTEGER :: row_start, row_end, col_start, col_end
  INTEGER :: II, JJ, KK, ilookup, jlookup
  INTEGER :: ind

  !! Get The Eigenvectors
  row_start = eigen_loop_start(1, exa%proc_rows, exa%rowid)
  row_end = eigen_loop_end(exa%mat_dim, exa%proc_rows, exa%rowid)
  col_start = eigen_loop_start(1, exa%proc_cols, exa%colid)
  col_end = eigen_loop_end(exa%mat_dim, exa%proc_cols, exa%colid)

  !! Convert to a 1D array for index ease.
  ALLOCATE(VD1(SIZE(VD, DIM = 1)*SIZE(VD, DIM = 2)))
  VD1(:) = PACK(VD, .TRUE.)

  !! First loop: count how many elements are needed
  !! originally I used Append to eliminate this loop but we need to save
  !! memory here.
  KK = 0
  ind = 1
  DO JJ = col_start, col_end
     DO II = row_start, row_end
        IF (ABS(VD1(ind + II -1)) .GT. params%threshold) THEN
           KK = KK + 1
        END IF
     END DO
     ind = ind + exa%offset
  END DO

  !! Construct The Triplet List
  CALL ConstructTripletList(triplet_v, KK)

  !! Reset indices and do the actual filling
  ind = 1
  KK = 0
  DO JJ = col_start, col_end
     jlookup = eigen_translate_l2g(JJ, exa%proc_cols, exa%colid)
     DO II = row_start, row_end
        IF (ABS(VD1(ind + II -1)) .GT. params%threshold) THEN
           KK = KK + 1
           ilookup = eigen_translate_l2g(II, exa%proc_rows, exa%rowid)
           CALL SetTriplet(triplet_v%DATA(KK), jlookup, ilookup, &
                & VD1(ind + II -1))
        END IF
     END DO
     ind = ind + exa%offset
  END DO
  DEALLOCATE(VD1)

  !! Fill and Clean Up
  CALL FillMatrixFromTripletList(V, triplet_v)
  CALL DestructTripletList(triplet_v)
