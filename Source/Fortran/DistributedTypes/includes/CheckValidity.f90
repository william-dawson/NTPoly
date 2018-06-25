  isvalid = .TRUE.
  !! Check allocation
  IF (.NOT. ALLOCATED(this%grid)) isvalid = .FALSE.
