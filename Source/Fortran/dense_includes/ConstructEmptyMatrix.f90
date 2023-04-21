  CALL DestructMatrix(this)

  this%rows = rows
  this%columns = columns

  ALLOCATE(this%DATA(rows, columns))
