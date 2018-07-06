  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, LMAT)
  ALLOCATE(local_norm(LMAT%columns))

  !! Sum Along Columns
  CALL MatrixColumnNorm(LMAT,local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
       & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

  !! Find Max Value Amonst Columns
  norm_value = MAXVAL(local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructMatrix(LMAT)
  DEALLOCATE(local_norm)
