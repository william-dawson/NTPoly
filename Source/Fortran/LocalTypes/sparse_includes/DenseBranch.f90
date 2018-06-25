  !! Handle Transposed Case
  IF (IsATransposed) THEN
     CALL untransposedMatA%Transpose(matA)
  ELSE
     CALL untransposedMatA%Copy(matA)
  END IF
  IF (IsBTransposed) THEN
     CALL untransposedMatB%Transpose(matB)
  ELSE
     CALL untransposedMatB%Copy(matB)
  END IF

  !! Convert Forward
  CALL ConvertSMatrixToD(untransposedMatA, DenseA)
  CALL ConvertSMatrixToD(untransposedMatB, DenseB)

  !! Multiply
  CALL MatrixMultiply(DenseA, DenseB, DenseC)

  !! Convert Back
  CALL ConvertDMatrixToS(DenseC, matC, threshold)
  CALL ScaleMatrix(matC,alpha)


  !! Cleanup
  CALL DenseA%Destruct
  CALL DenseB%Destruct
  CALL DenseC%Destruct
