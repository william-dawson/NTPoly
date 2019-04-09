  CALL ConstructEmptyMatrix(matAT, matA%columns, matA%rows)
  matAT%data = TRANSPOSE(matA%data)
