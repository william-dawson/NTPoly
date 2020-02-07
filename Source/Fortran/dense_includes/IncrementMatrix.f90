  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0_NTREAL
  ELSE
     alpha = alpha_in
  END IF

  MatB%data = MatB%data + alpha*MatA%data
