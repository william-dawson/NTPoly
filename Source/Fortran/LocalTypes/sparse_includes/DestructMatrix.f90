  IF (ALLOCATED(this%outer_index)) THEN
     DEALLOCATE(this%outer_index)
  END IF
  IF (ALLOCATED(this%inner_index)) THEN
     DEALLOCATE(this%inner_index)
  END IF
  IF (ALLOCATED(this%values)) THEN
     DEALLOCATE(this%values)
  END IF
