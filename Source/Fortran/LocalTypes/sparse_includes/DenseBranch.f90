  !! Handle Transposed Case
  IF (IsATransposed) THEN
     CALL TransposeMatrix(matA,untransposedMatA)
  ELSE
     CALL CopyMatrix(matA,untransposedMatA)
  END IF
  IF (IsBTransposed) THEN
     CALL TransposeMatrix(matB,untransposedMatB)
  ELSE
     CALL CopyMatrix(matB,untransposedMatB)
  END IF

  !! Convert Forward
  CALL ConstructMatrixDFromS(untransposedMatA, DenseA)
  CALL ConstructMatrixDFromS(untransposedMatB, DenseB)

  !! Multiply
  CALL MultiplyMatrix(DenseA, DenseB, DenseC)

  !! Convert Back
  CALL ConstructMatrixSFromD(DenseC, matC, threshold)
  CALL ScaleMatrix(matC,alpha)

  !! Cleanup
  CALL DestructMatrix(DenseA)
  CALL DestructMatrix(DenseB)
  CALL DestructMatrix(DenseC)
