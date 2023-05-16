  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, LMAT)

  ALLOCATE(column_sigma_contribution(LMAT%columns))
  column_sigma_contribution = 0
  DO II = 1, LMAT%columns
     DO JJ = LMAT%outer_index(II), LMAT%outer_index(II + 1) - 1
        column_sigma_contribution(II) = column_sigma_contribution(II) + &
             & ABS(LMAT%values(JJ + 1))
     END DO
  END DO
  CALL MPI_Allreduce(MPI_IN_PLACE, column_sigma_contribution, &
       & LMAT%columns, MPINTREAL, MPI_SUM, &
       & this%process_grid%column_comm, ierr)
  CALL MPI_Allreduce(MAXVAL(column_sigma_contribution), sigma_value, 1, &
       & MPINTREAL, MPI_MAX, this%process_grid%row_comm, ierr)
  sigma_value = 1.0_NTREAL / (sigma_value**2)

  DEALLOCATE(column_sigma_contribution)
  CALL DestructMatrix(LMAT)
