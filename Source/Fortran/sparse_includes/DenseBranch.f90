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
  ! Match SparseBranch: threshold the scaled product, not the raw AB.
  DenseC%DATA = alpha*DenseC%DATA
  CALL ConstructMatrixSFromD(DenseC, matC, threshold)

  !! Cleanup
  CALL DestructMatrix(DenseC)
