  CALL this%Destruct
  ALLOCATE(this%outer_index, source=matA%outer_index)
  ALLOCATE(this%inner_index, source=matA%inner_index)
  ALLOCATE(this%values, source=matA%values)
  this%rows = matA%rows
  this%columns = matA%columns
