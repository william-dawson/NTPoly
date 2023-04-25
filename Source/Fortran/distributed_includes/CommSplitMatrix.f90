  !! For Grid Splitting
  TYPE(ProcessGrid_t) :: new_grid
  INTEGER :: between_grid_comm
  INTEGER :: between_grid_size
  INTEGER :: between_grid_rank
  !! For Data Redistribution
  INTEGER :: fsize
  INTEGER :: counter
  INTEGER :: ierr

  IF (this%process_grid%total_processors .EQ. 1) THEN
     CALL CopyMatrix(this, split_mat)
     my_color = 0
     split_slice = .TRUE.
  ELSE
     !! Split The Grid
     CALL SplitProcessGrid(this%process_grid, new_grid, my_color, &
          & split_slice, between_grid_comm)

     !! Copy The Data Across New Process Grids. Unnecessary if we just split
     !! by slices.
     CALL GetMatrixTripletList(this, full_list)
     IF (.NOT. split_slice) THEN
        CALL MPI_COMM_SIZE(between_grid_comm, between_grid_size, ierr)
        CALL MPI_COMM_RANK(between_grid_comm, between_grid_rank, ierr)

        !! Build Send Lists
        fsize = full_list%CurrentSize
        ALLOCATE(send_list(between_grid_size))
        IF (my_color .EQ. 0) THEN
           !! The smaller process grid only needs to send to process 2
           CALL ConstructTripletList(send_list(1))
           CALL ConstructTripletList(send_list(2), full_list%CurrentSize)
           send_list(2)%DATA(:fsize) = full_list%DATA(:fsize)
           DO counter = 3, between_grid_size
              CALL ConstructTripletList(send_list(counter))
           END DO
        ELSE
           !! The larger process grid only needs to send to process 1
           CALL ConstructTripletList(send_list(1), full_list%CurrentSize)
           send_list(1)%DATA(:fsize) = full_list%DATA(:fsize)
           DO counter = 2, between_grid_size
              CALL ConstructTripletList(send_list(counter))
           END DO
        END IF
        CALL ConstructTripletList(send_list(between_grid_rank + 1), &
             & full_list%CurrentSize)
        send_list(between_grid_rank + 1)%DATA(:fsize) = full_list%DATA(:fsize)
        CALL RedistributeTripletLists(send_list, between_grid_comm, new_list)
     END IF

     !! Create The New Matrix
     CALL ConstructEmptyMatrix(split_mat, this%actual_matrix_dimension, &
          & process_grid_in=new_grid, is_complex_in=this%is_complex)
     IF (.NOT. split_slice) THEN
        CALL FillMatrixFromTripletList(split_mat, new_list, .TRUE.)
     ELSE
        CALL FillMatrixFromTripletList(split_mat, full_list, .TRUE.)
     END IF

     !! Cleanup
     CALL DestructTripletList(full_list)
     CALL DestructTripletList(new_list)
     IF (ALLOCATED(send_list)) THEN
        DO counter = 1, between_grid_size
           CALL DestructTripletList(send_list(counter))
        END DO
     END IF
  END IF
