  CALL matB%Destruct
  ALLOCATE(matB%outer_index, source=matA%outer_index)
  ALLOCATE(matB%inner_index, source=matA%inner_index)
  ALLOCATE(matB%values, source=matA%values)
  matB%rows = matA%rows
  matB%columns = matA%columns
