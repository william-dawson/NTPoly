  CALL ConstructEmptyMatrix(matAT, matA%rows, matA%columns)
  matAT%data = TRANSPOSE(matA%data)
