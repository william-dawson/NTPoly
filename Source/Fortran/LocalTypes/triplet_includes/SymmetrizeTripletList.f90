  !! Local variables
  INTEGER :: counter
  INTEGER :: initial_size

  initial_size = this%CurrentSize
  DO counter = 1, initial_size
     CALL this%Get(counter,temporary)
     IF (temporary%index_row .EQ. temporary%index_column) CYCLE

     CALL temporary%Transpose()
     
     SELECT CASE(pattern_type)
     CASE(MM_HERMITIAN)
        CALL temporary%Conjg()
        CALL this%Append(temporary)
     CASE(MM_SKEW_SYMMETRIC)
        CALL temporary%Scale(REAL(-1.0,NTREAL))
        CALL this%Append(temporary)
      CASE(MM_SYMMETRIC)
        CALL this%Append(temporary)
     END SELECT
  END DO
