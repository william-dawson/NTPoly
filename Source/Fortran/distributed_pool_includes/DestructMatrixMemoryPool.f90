  !! Allocate
  IF (ALLOCATED(this%grid)) THEN
     DO II = LBOUND(this%grid, 2), UBOUND(this%grid, 2)
        DO JJ = LBOUND(this%grid, 1), UBOUND(this%grid, 1)
           CALL DestructMatrixMemoryPool(this%grid(JJ, II))
        END DO
     END DO
     DEALLOCATE(this%grid)
  END IF
