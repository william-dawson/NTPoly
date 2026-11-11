! Independent expected values across the dense/sparse local GEMM crossover.
PROGRAM TestGemmThreshold
  USE DataTypesModule, ONLY: NTREAL
  USE SMatrixModule
  USE SMatrixAlgebraModule
  USE TripletModule
  USE TripletListModule
  IMPLICIT NONE
  TYPE(Matrix_lsr) :: a,b,c
  TYPE(TripletList_r) :: entries
  TYPE(Triplet_r) :: entry
  INTEGER :: n,j,trial
  REAL(NTREAL) :: alpha,threshold,expected
  DO n=4,20,16
    CALL ConstructTripletList(entries)
    DO j=1,n
      entry%index_row=j; entry%index_column=j; entry%point_value=.25_NTREAL
      CALL AppendToTripletList(entries,entry)
    END DO
    CALL ConstructMatrixFromTripletList(a,entries,n,n)
    CALL DestructTripletList(entries); CALL CopyMatrix(a,b)
    CALL ScaleMatrix(b,4.0_NTREAL)
    DO trial=1,5
      SELECT CASE(trial)
      CASE(1); alpha=2; threshold=.3_NTREAL
      CASE(2); alpha=.5_NTREAL; threshold=.2_NTREAL
      CASE(3); alpha=0; threshold=0
      CASE(4); alpha=-2; threshold=.3_NTREAL
      CASE(5); alpha=2; threshold=.5_NTREAL
      END SELECT
      CALL MatrixMultiply(a,b,c,alpha_in=alpha,threshold_in=threshold)
      expected=alpha*.25_NTREAL
      IF(ABS(expected)>threshold) THEN
        IF(SIZE(c%values)/=n) ERROR STOP 'missing scaled terms'
        IF(ANY(c%values/=expected)) ERROR STOP 'incorrect scaled values'
      ELSE
        IF(SIZE(c%values)/=0) ERROR STOP 'unpruned scaled terms'
      END IF
      CALL DestructMatrix(c)
    END DO
    CALL DestructMatrix(a); CALL DestructMatrix(b)
  END DO
  PRINT *, 'Dense/sparse GEMM scaled-threshold consistency: PASS'
END PROGRAM
