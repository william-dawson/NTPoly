  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL IncrementMatrix(matA%LMAT(II, JJ), matB%LMAT(II, JJ), alpha, &
             & threshold)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
