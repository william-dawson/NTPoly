  isvalid = .TRUE.
  !! Check allocation
  IF (.NOT. ALLOCATED(this%pruned_list)) isvalid = .FALSE.
  IF (.NOT. ALLOCATED(this%value_array)) isvalid = .FALSE.

  IF (isvalid) THEN
     !! Check allocation size
     IF (.NOT. SIZE(this%value_array, dim=2) .EQ. rows) THEN
        isvalid = .FALSE.
     END IF
     IF (.NOT. SIZE(this%value_array, dim=1) .EQ. columns) THEN
        isvalid = .FALSE.
     END IF
  END IF
