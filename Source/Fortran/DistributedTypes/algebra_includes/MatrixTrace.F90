  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, LMAT)

  !! Compute The Local Contribution
  trace_value = 0
  CALL MatrixToTripletList(LMAT, TLIST)
  DO counter = 1, TLIST%CurrentSize
     IF (this%start_row + TLIST%data(counter)%index_row .EQ. &
          & this%start_column + TLIST%data(counter)%index_column) THEN
        trace_value = trace_value + &
             & REAL(TLIST%data(counter)%point_value, NTREAL)
     END IF
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)

  CALL DestructMatrix(LMAT)
