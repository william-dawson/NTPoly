  CALL ConstructEmptyMatrix(matC, matA%actual_matrix_dimension, &
       & matA%process_grid, matA%is_complex)

  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL PairwiseMultiplyMatrix(matA%LMAT(II, JJ), matB%LMAT(II, JJ), &
             & matC%LMAT(II, JJ))
     END DO
  END DO
  !$omp end do
  !$omp end parallel
