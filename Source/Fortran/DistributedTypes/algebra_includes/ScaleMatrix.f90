  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL ScaleMatrix(this%LOCALDATA(II,JJ),constant)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
