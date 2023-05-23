  CALL ConstructEmptyMatrix(matB, matA%rows, matA%columns)
  matB%DATA(:, :) = matA%DATA
