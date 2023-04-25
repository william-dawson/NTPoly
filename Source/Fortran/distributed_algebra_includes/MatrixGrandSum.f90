  sum = 0
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL MatrixGrandSum(this%LMAT(II, JJ), TEMP)
        sum = sum + REAL(TEMP, KIND=NTREAL)
     END DO
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPIDATATYPE, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
