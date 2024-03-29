  !! Allocate Space For Result
  ALLOCATE(per_column_min(this%local_columns))
  ALLOCATE(per_column_max(this%local_columns))

  !! Compute The Local Contribution
  per_column_min = 0
  per_column_max = 0
  CALL GetMatrixTripletList(this, tlist)
  DO II = 1, tlist%CurrentSize
     local_column = tlist%DATA(II)%index_column - &
          & this%start_column + 1
     IF (tlist%DATA(II)%index_row .EQ. tlist%DATA(II)%index_column) THEN
        per_column_min(local_column) = per_column_min(local_column) + &
             & REAL(tlist%DATA(II)%point_value, KIND = NTREAL)
        per_column_max(local_column) = per_column_max(local_column) + &
             & REAL(tlist%DATA(II)%point_value, KIND = NTREAL)
     ELSE
        per_column_min(local_column) = per_column_min(local_column) - &
             & ABS(tlist%DATA(II)%point_value)
        per_column_max(local_column) = per_column_max(local_column) + &
             & ABS(tlist%DATA(II)%point_value)
     END IF
  END DO

  !! Sum Along Columns
  CALL MPI_Allreduce(MPI_IN_PLACE, per_column_min, SIZE(per_column_min), &
       & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE, per_column_max, SIZE(per_column_max), &
       & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

  min_value = MINVAL(per_column_min)
  max_value = MAXVAL(per_column_max)

  CALL MPI_Allreduce(MPI_IN_PLACE, min_value, 1, MPINTREAL, MPI_MIN, &
       & this%process_grid%row_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE, max_value, 1, MPINTREAL, MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructTripletList(tlist)
  DEALLOCATE(per_column_min)
  DEALLOCATE(per_column_max)
