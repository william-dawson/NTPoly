  INTEGER :: rows_a, rows_b, rows_c
  INTEGER :: cols_a, cols_b, cols_c
  INTEGER :: II, JJ, KK

  !! Get the blocking parameters
  rows_a = GetMatrixBlockRows(matA)
  rows_b = GetMatrixBlockRows(matB)
  rows_c = GetMatrixBlockRows(matC)