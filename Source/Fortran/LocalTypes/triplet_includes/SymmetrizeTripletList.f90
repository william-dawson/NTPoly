  !! Local variables
  INTEGER :: counter
  INTEGER :: initial_size

  initial_size = this%CurrentSize
  DO counter = 1, initial_size
     CALL this%Get(counter,temporary)
     IF (temporary%index_row .EQ. temporary%index_column) CONTINUE
     CALL temporary%Transpose()

     SELECT CASE(pattern_type)
     CASE(MM_HERMITIAN)
        CALL temporary%Conjg()
     CASE(MM_SKEW_SYMMETRIC)
        CALL temporary%Scale(REAL(-1.0,NTREAL))
     END SELECT

     CALL this%Append(temporary)
  END DO
