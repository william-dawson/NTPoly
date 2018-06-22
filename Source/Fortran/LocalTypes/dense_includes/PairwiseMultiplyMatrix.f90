  !! Local Variables
  INTEGER :: II, JJ

  CALL matC%InitEmpty(matA%rows, matA%columns)
  DO JJ = 1, matA%columns
    DO II = 1, matA%rows
       matC%data(II,JJ) = matA%data(II,JJ)*matB%data(II,JJ)
    END DO
  END DO
