  IF (this%process_grid%number_of_blocks_rows == 1 .AND. &
       & this%process_grid%number_of_blocks_columns == 1) THEN
     ! Intrinsic assignment owns an independent copy of the allocatable arrays.
     ! No split, temporary transpose or reordering is needed for one block.
     LOCALMATRIX(1,1) = matrix_to_split
     ! Empty constructors may omit these arrays; block storage requires them.
     IF (.NOT. ALLOCATED(LOCALMATRIX(1,1)%inner_index)) &
          & ALLOCATE(LOCALMATRIX(1,1)%inner_index(0))
     IF (.NOT. ALLOCATED(LOCALMATRIX(1,1)%values)) &
          & ALLOCATE(LOCALMATRIX(1,1)%values(0))
  ELSE
     CALL SplitMatrix(matrix_to_split, &
          & this%process_grid%number_of_blocks_rows, &
          & this%process_grid%number_of_blocks_columns, LOCALMATRIX)
  END IF
