  INTEGER :: II, JJ

  IF (ALLOCATED(this%h_data)) THEN
     DO II = 1, SIZE(this%h_data, DIM=1)
     	DO JJ = 1, SIZE(this%h_data, DIM=2)
     		CALL DestructMatrix(this%h_data(II,JJ))
     	END DO
     END DO 
  END IF

  IF (ALLOCATED(this%base_data)) THEN
     DO II = 1, SIZE(this%base_data, DIM=1)
     	DO JJ = 1, SIZE(this%base_data, DIM=2)
     		CALL DestructMatrix(this%base_data(II,JJ))
     	END DO
     END DO 
  END IF