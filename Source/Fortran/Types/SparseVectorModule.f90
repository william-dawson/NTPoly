!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
!! This module can add two of those vectors together.
MODULE SparseVectorModule
  USE DataTypesModule, ONLY : NTREAL
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AddSparseVectors
  PUBLIC :: DotSparseVectors
  PUBLIC :: PairwiseMultiplyVectors
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add together two sparse vectors. C = A + alpha*B
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  !! @param[out] inner_index_c list of indices computed for C.
  !! @param[out] values_c list of values computed for C.
  !! @param[out] total_values_c this is the total number of values in C.
  !! @param[in] alpha_in value to scale VecB by. Optional, default is 1.0.
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! The values that are returned for C are only valid in the range
  !! (1:total_values_c). We do not do an automatic shrinking of the array
  !! to keep this routine low in overhead.
  !! @todo in principal this can be done without any branching.
  PURE SUBROUTINE AddSparseVectors(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_a
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_b
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Data
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    !! Temporary Variables
    INTEGER :: working_index_a, working_index_b
    REAL(NTREAL) :: working_value_a, working_value_b
    !! Counter Variables
    INTEGER :: counter_a, counter_b, counter_c

    !! Process Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0d+0
    ELSE
       threshold = threshold_in
    END IF

    counter_a = 1
    counter_b = 1
    counter_c = 1
    sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
         & SIZE(inner_index_b))
       !! Current inner indices and values
       working_index_a = inner_index_a(counter_a)
       working_index_b = inner_index_b(counter_b)
       working_value_a = alpha*values_a(counter_a)
       working_value_b = values_b(counter_b)
       !! Figure out from which vector an insertion will be performed
       IF (working_index_a .EQ. working_index_b) THEN
          IF (ABS(working_value_a + working_value_b) .GT. threshold) THEN
             inner_index_c(counter_c) = working_index_a
             values_c(counter_c) = working_value_a + working_value_b
             counter_c = counter_c + 1
          END IF
          counter_a = counter_a + 1
          counter_b = counter_b + 1
       ELSE IF (working_index_a .GT. working_index_b) THEN
          IF (ABS(working_value_b) .GT. threshold) THEN
             inner_index_c(counter_c) = working_index_b
             values_c(counter_c) = working_value_b
             counter_c = counter_c + 1
          END IF
          counter_b = counter_b + 1
       ELSE !! implies working_index_b > working_index_b
          IF (ABS(working_value_a) .GT. threshold) THEN
             inner_index_c(counter_c) = working_index_a
             values_c(counter_c) = working_value_a
             counter_c = counter_c + 1
          END IF
          counter_a = counter_a + 1
       END IF
    END DO sum_loop

    !! Handle case where one was blank
    cleanup_a: DO WHILE (counter_a .LE. SIZE(inner_index_a))
       inner_index_c(counter_c) = inner_index_a(counter_a)
       values_c(counter_c) = values_a(counter_a)*alpha
       counter_a = counter_a + 1
       counter_c = counter_c + 1
    END DO cleanup_a
    cleanup_b: DO WHILE (counter_b .LE. SIZE(inner_index_b))
       inner_index_c(counter_c) = inner_index_b(counter_b)
       values_c(counter_c) = values_b(counter_b)
       counter_b = counter_b + 1
       counter_c = counter_c + 1
    END DO cleanup_b

    total_values_c = counter_c - 1
  END SUBROUTINE AddSparseVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  PURE FUNCTION DotSparseVectors(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_a
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_b
    REAL(NTREAL) :: product
    !! Temporary Variables
    INTEGER :: working_index_a, working_index_b
    REAL(NTREAL) :: working_value_a, working_value_b
    !! Counter Variables
    INTEGER :: counter_a, counter_b

    counter_a = 1
    counter_b = 1
    product = 0
    sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
         & SIZE(inner_index_b))
       !! Current inner indices and values
       working_index_a = inner_index_a(counter_a)
       working_index_b = inner_index_b(counter_b)
       working_value_a = values_a(counter_a)
       working_value_b = values_b(counter_b)
       !! Figure out from which vector an insertion will be performed
       IF (working_index_a .EQ. working_index_b) THEN
          product = product + working_value_a * working_value_b
          counter_a = counter_a + 1
          counter_b = counter_b + 1
       ELSE IF (working_index_a .GT. working_index_b) THEN
          counter_b = counter_b + 1
       ELSE !! implies working_index_b > working_index_b
          counter_a = counter_a + 1
       END IF
    END DO sum_loop

  END FUNCTION DotSparseVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  !! @param[out] inner_index_c list of indices computed for C.
  !! @param[out] values_c list of values computed for C.
  !! @param[out] total_values_c this is the total number of values in C.
  PURE SUBROUTINE PairwiseMultiplyVectors(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_a
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_b
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    INTEGER :: working_index_a, working_index_b
    REAL(NTREAL) :: working_value_a, working_value_b
    !! Counter Variables
    INTEGER :: counter_a, counter_b, counter_c

    counter_a = 1
    counter_b = 1
    counter_c = 1
    sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
         & SIZE(inner_index_b))
       !! Current inner indices and values
       working_index_a = inner_index_a(counter_a)
       working_index_b = inner_index_b(counter_b)
       working_value_a = values_a(counter_a)
       working_value_b = values_b(counter_b)
       !! Figure out from which vector an insertion will be performed
       IF (working_index_a .EQ. working_index_b) THEN
          inner_index_c(counter_c) = working_index_a
          values_c(counter_c) = working_value_a * working_value_b
          counter_c = counter_c + 1
          counter_a = counter_a + 1
          counter_b = counter_b + 1
       ELSE IF (working_index_a .GT. working_index_b) THEN
          counter_b = counter_b + 1
       ELSE !! implies working_index_b > working_index_b
          counter_a = counter_a + 1
       END IF
    END DO sum_loop
    total_values_c = counter_c - 1
  END SUBROUTINE PairwiseMultiplyVectors
END MODULE SparseVectorModule
