  !! Temporary copy
  ALLOCATE(temporary_data(this%CurrentSize))
  temporary_data = this%DATA(:this%CurrentSize)

  !! Create new memory
  IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
  ALLOCATE(this%DATA(size))

  !! Copy back
  this%DATA(:this%CurrentSize) = temporary_data

  !! Cleanup
  DEALLOCATE(temporary_data)
