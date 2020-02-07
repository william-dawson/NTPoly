  IF (ALLOCATED(this%h_data)) THEN
     base_case = .FALSE.
  ELSE
     base_case = .TRUE.
  END IF