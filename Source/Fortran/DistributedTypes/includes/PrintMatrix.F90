  !! Helpers For Communication
  TYPE(ReduceHelper_t) :: row_helper
  TYPE(ReduceHelper_t) :: column_helper
  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: ierr

  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data)

  !! Merge Columns
  CALL TransposeMatrix(merged_local_data, merged_local_dataT)
  CALL ReduceMatrixSizes(merged_local_dataT, this%process_grid%column_comm, &
       & column_helper)
  CALL MPI_Wait(column_helper%size_request ,mpi_status,ierr)
  CALL ReduceAndComposeMatrixData(merged_local_dataT, &
       & this%process_grid%column_comm,merged_columns, &
       & column_helper)
  CALL MPI_Wait(column_helper%outer_request,mpi_status,ierr)
  CALL MPI_Wait(column_helper%inner_request,mpi_status,ierr)
  CALL MPI_Wait(column_helper%data_request,mpi_status,ierr)
  CALL ReduceAndComposeMatrixCleanup(merged_local_dataT,merged_columns, &
       & column_helper)

  !! Merge Rows
  CALL TransposeMatrix(merged_columns,merged_columnsT)
  CALL ReduceMatrixSizes(merged_columnsT, this%process_grid%row_comm, row_helper)
  CALL MPI_Wait(row_helper%size_request,mpi_status,ierr)
  CALL ReduceAndComposeMatrixData(merged_columnsT, this%process_grid%row_comm, &
       & full_gathered,row_helper)
  CALL MPI_Wait(row_helper%outer_request,mpi_status,ierr)
  CALL MPI_Wait(row_helper%inner_request,mpi_status,ierr)
  CALL MPI_Wait(row_helper%data_request,mpi_status,ierr)
  CALL ReduceAndComposeMatrixCleanup(merged_columnsT,full_gathered,row_helper)

  !! Make these changes so that it prints the logical rows/columns
  full_gathered%rows = this%actual_matrix_dimension
  full_gathered%columns = this%actual_matrix_dimension

  IF (IsRoot(this%process_grid)) THEN
     IF (PRESENT(file_name_in)) THEN
        CALL PrintMatrix(full_gathered, file_name_in)
     ELSE
        CALL PrintMatrix(full_gathered)
     END IF
  END IF

  CALL DestructMatrix(merged_local_data)
  CALL DestructMatrix(merged_local_dataT)
  CALL DestructMatrix(merged_columns)
  CALL DestructMatrix(merged_columnsT)
