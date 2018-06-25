  !! Allocate
  IF (ALLOCATED(this%grid)) THEN
     DO column_counter = LBOUND(this%grid,2), UBOUND(this%grid,2)
        DO row_counter = LBOUND(this%grid,1), UBOUND(this%grid,1)
           CALL this%grid(row_counter, column_counter)%Destruct
        END DO
     END DO
     DEALLOCATE(this%grid)
  END IF
