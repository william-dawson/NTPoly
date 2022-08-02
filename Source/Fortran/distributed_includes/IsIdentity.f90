  !! Local Data
  INTEGER :: II
  INTEGER :: ierr
  INTEGER :: diag_count

  is_identity = .TRUE.
  diag_count = 0
  CALL GetMatrixTripletList(this, tlist)
  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, trip)
     IF (trip%index_row .NE. trip%index_column) THEN
        is_identity = .FALSE.
        EXIT
     ELSE IF (trip%point_value .NE. 1.0_NTREAL) THEN
        is_identity = .FALSE.
     ELSE
        diag_count = diag_count + 1
     END IF
  END DO

  CALL MPI_Allreduce(MPI_IN_PLACE, is_identity, 1, MPI_LOGICAL, &
       & MPI_LAND, this%process_grid%within_slice_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE, diag_count, 1, MPINTINTEGER, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
  IF (.NOT. is_identity .OR. &
       & .NOT. (diag_count .EQ. this%logical_matrix_dimension)) THEN
     is_identity = .FALSE.
  END IF
