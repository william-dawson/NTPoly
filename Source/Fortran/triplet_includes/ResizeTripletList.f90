  INTEGER :: old_size

  !! Temporary copy
  old_size = this%CurrentSize
  ALLOCATE(temporary_data(old_size))
  temporary_data(:) = this%DATA(:old_size)

  !! Create new memory
  IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
  ALLOCATE(this%DATA(size))

  !! Copy back
  IF (old_size .LT. size) THEN
     this%DATA(:old_size) = temporary_data(:old_size)
  ELSE
     this%DATA(:size) = temporary_data(:size)
  END IF

  !! Cleanup
  DEALLOCATE(temporary_data)
