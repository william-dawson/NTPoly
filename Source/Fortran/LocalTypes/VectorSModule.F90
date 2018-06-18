!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
MODULE VectorSModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AddSparseVectors
  PUBLIC :: DotSparseVectors
  PUBLIC :: PairwiseMultiplyVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE AddSparseVectors
     MODULE PROCEDURE AddSparseVectors_r
     MODULE PROCEDURE AddSparseVectors_c
  END INTERFACE
  INTERFACE DotSparseVectors
     MODULE PROCEDURE DotSparseVectors_r
     MODULE PROCEDURE DotSparseVectors_c
  END INTERFACE
  INTERFACE PairwiseMultiplyVectors
     MODULE PROCEDURE PairwiseMultiplyVectors_r
     MODULE PROCEDURE PairwiseMultiplyVectors_c
  END INTERFACE
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
  PURE SUBROUTINE AddSparseVectors_r(inner_index_a,values_a,inner_index_b, &
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
    !! Local Variables
    REAL(NTREAL) :: working_value_a, working_value_b

    include "sparse_includes/AddSparseVectors.f90"

  END SUBROUTINE AddSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  PURE FUNCTION DotSparseVectors_r(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_a
    REAL(NTREAL), DIMENSION(:), INTENT(IN)    :: values_b
    REAL(NTREAL) :: product
    !! Temporary Variables
    REAL(NTREAL) :: working_value_a, working_value_b

    include "sparse_includes/DotSparseVectors.f90"

  END FUNCTION DotSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  !! @param[out] inner_index_c list of indices computed for C.
  !! @param[out] values_c list of values computed for C.
  !! @param[out] total_values_c this is the total number of values in C.
  PURE SUBROUTINE PairwiseMultiplyVectors_r(inner_index_a,values_a,inner_index_b,&
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
    REAL(NTREAL) :: working_value_a, working_value_b

    include "sparse_includes/PairwiseMultiplyVectors.f90"

  END SUBROUTINE PairwiseMultiplyVectors_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  PURE SUBROUTINE AddSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_a
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_b
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    include "sparse_includes/AddSparseVectors.f90"

  END SUBROUTINE AddSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  PURE FUNCTION DotSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_a
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_b
    COMPLEX(NTCOMPLEX) :: product
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    include "sparse_includes/DotSparseVectors.f90"

  END FUNCTION DotSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  !! @param[in] inner_index_a list of indices for A.
  !! @param[in] values_a list of values for A.
  !! @param[in] inner_index_b list of indices for B.
  !! @param[in] values_b list of values for B.
  !! @param[out] inner_index_c list of indices computed for C.
  !! @param[out] values_c list of values computed for C.
  !! @param[out] total_values_c this is the total number of values in C.
  PURE SUBROUTINE PairwiseMultiplyVectors_c(inner_index_a,values_a,inner_index_b,&
       & values_b,inner_index_c,values_c,total_values_c)
    !! Parameters
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_a
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)    :: values_b
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    include "sparse_includes/PairwiseMultiplyVectors.f90"

  END SUBROUTINE PairwiseMultiplyVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE VectorSModule
