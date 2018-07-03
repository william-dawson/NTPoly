  !! Local Data
  TYPE(Matrix_ps) :: matAH
  TYPE(Matrix_ps) :: matC

  IF (matA%is_complex) THEN
     CALL CopyMatrix(matA, matAH)
     CALL ConjugateMatrix(matAH)
     CALL PairwiseMultiplyMatrix(matAH, matB, matC)
     CALL DestructMatrix(matAH)
  ELSE
     CALL PairwiseMultiplyMatrix(matA,matB,matC)
  END IF

  CALL MatrixGrandSum(matC, product)
  CALL DestructMatrix(matC)
