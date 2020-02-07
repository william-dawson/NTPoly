  CALL ComposeMatrix(this, merged)
  CALL MatrixColumnNorm(merged, norm_per_column)
  CALL DestructMatrix(merged)
