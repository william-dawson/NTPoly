  !! Local Data
  INTEGER :: II
  INTEGER :: ierr
  INTEGER :: diag_count

  !! Loop over elements checking that they are on the diagonal and are 1
  is_identity = .TRUE.
  diag_count = 0
  CALL GetMatrixTripletList(this, tlist)
  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, trip)
     IF (trip%index_row .NE. trip%index_column) THEN
        is_identity = .FALSE.
        EXIT
#ifdef ISCOMPLEX
     ELSE IF (ABS(trip%point_value - 1.0_NTCOMPLEX) &
              & .GT. TINY(trip%point_value)) THEN
#else
     ELSE IF (REAL(trip%point_value - 1.0_NTREAL) .GT. TINY(1.0_NTREAL) .OR. &
            & AIMAG(trip%point_value - 1.0_NTREAL) .GT. TINY(1.0_NTREAL)) THEN
#endif
        is_identity = .FALSE.
     ELSE
        diag_count = diag_count + 1
     END IF
  END DO

  !! Share findings across processes
  CALL MPI_Allreduce(MPI_IN_PLACE, is_identity, 1, MPI_LOGICAL, &
       & MPI_LAND, this%process_grid%within_slice_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE, diag_count, 1, MPINTINTEGER, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)

  !! Make sure we found enough 1s
  IF (.NOT. is_identity .OR. &
       & .NOT. (diag_count .EQ. this%actual_matrix_dimension)) THEN
     is_identity = .FALSE.
  END IF
