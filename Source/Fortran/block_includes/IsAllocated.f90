  is_allocated = .FALSE.
  IF (ALLOCATED(this%h_data)) THEN
     is_allocated = .TRUE.
  ELSEIF (ALLOCATED(this%base_data)) THEN
     is_allocated = .TRUE.
  END IF