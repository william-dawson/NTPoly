  !! Temporary copy
  ALLOCATE(temporary_data(this%CurrentSize))
  temporary_data = this%data(:this%CurrentSize)

  !! Create new memory
  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  ALLOCATE(this%data(size))

  !! Copy back
  this%data(:this%CurrentSize) = temporary_data

  !! Cleanup
  DEALLOCATE(temporary_data)
