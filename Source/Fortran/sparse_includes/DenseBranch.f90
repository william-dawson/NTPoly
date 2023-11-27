  !! Convert Forward
  CALL ConstructMatrixDFromS(matA, DenseA)
  CALL ConstructMatrixDFromS(matB, DenseB)

  !! Multiply
  CALL MultiplyMatrix(DenseA, DenseB, DenseC, &
       & IsATransposed_in = IsATransposed, IsBTransposed_in = IsBTransposed)

  !! Cleanup Intermediate
  CALL DestructMatrix(DenseA)
  CALL DestructMatrix(DenseB)

  !! Convert Back
  CALL ConstructMatrixSFromD(DenseC, matC, threshold)
  CALL ScaleMatrix(matC, alpha)

  !! Cleanup
  CALL DestructMatrix(DenseC)
