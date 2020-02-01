  CALL ConstructMatrixDFromS(matA, DenseA)
  CALL ConstructMatrixDFromS(matB, DenseB)

  !! Handle Transposed Case
  IF (IsATransposed) THEN
     CALL TransposeMatrix(DenseA,untransposedMatA)
  ELSE
     CALL CopyMatrix(DenseA,untransposedMatA)
  END IF
  IF (IsBTransposed) THEN
     CALL TransposeMatrix(DenseB,untransposedMatB)
  ELSE
     CALL CopyMatrix(DenseB,untransposedMatB)
  END IF

  !! Convert Forward

  !! Multiply
  CALL MultiplyMatrix(untransposedMatA, untransposedMatB, DenseC)

  !! Convert Back
  CALL ConstructMatrixSFromD(DenseC, matC, threshold)
  CALL ScaleMatrix(matC,alpha)

  !! Cleanup
  CALL DestructMatrix(untransposedMatA)
  CALL DestructMatrix(untransposedMatB)
  CALL DestructMatrix(DenseA)
  CALL DestructMatrix(DenseB)
  CALL DestructMatrix(DenseC)
