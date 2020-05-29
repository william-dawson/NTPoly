  CALL ConstructEmptyMatrix(matAT, matA%columns, matA%rows)
  matAT%DATA = TRANSPOSE(matA%DATA)
