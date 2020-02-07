  REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: temp
  INTEGER :: cols

  cols = GetMatrixColumns(this)
  ALLOCATE(temp(cols))
  norm = MAXVAL(temp)
  DEALLOCATE(temp)