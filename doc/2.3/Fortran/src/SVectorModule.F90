!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!> Compressed vectors are stored in two lists. The first is a list of indices,
!> the second a list of values.
MODULE SVectorModule
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
  !> The values that are returned for C are only valid in the range
  !> (1:total_values_c). We do not do an automatic shrinking of the array
  !> to keep this routine low in overhead.
  PURE SUBROUTINE AddSparseVectors_r(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    !> Value to scale VecB by. Optional, default is 1.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> for flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> The total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Variables
    REAL(NTREAL) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/AddSparseVectors.f90"

  END SUBROUTINE AddSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add together two sparse vectors. C = A + alpha*B
  !> The values that are returned for C are only valid in the range
  !> (1:total_values_c). We do not do an automatic shrinking of the array
  !> to keep this routine low in overhead.
  PURE SUBROUTINE AddSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    !> Value to scale VecB by. Optional, default is 1.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> for flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> The total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/AddSparseVectors.f90"

  END SUBROUTINE AddSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  PURE FUNCTION DotSparseVectors_r(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_b
    !> Dot product.
    REAL(NTREAL) :: product
    !! Temporary Variables
    REAL(NTREAL) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/DotSparseVectors.f90"

  END FUNCTION DotSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  PURE FUNCTION DotSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN) :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN) :: values_b
    !> Dot product.
    COMPLEX(NTCOMPLEX) :: product
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/DotSparseVectors.f90"

  END FUNCTION DotSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  PURE SUBROUTINE PairwiseMultiplyVectors_r(inner_index_a,values_a, &
       & inner_index_b,values_b,inner_index_c,values_c,total_values_c)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices computed for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    !> This is the total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    REAL(NTREAL) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/PairwiseMultiplyVectors.f90"

  END SUBROUTINE PairwiseMultiplyVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  PURE SUBROUTINE PairwiseMultiplyVectors_c(inner_index_a,values_a, &
       & inner_index_b, values_b,inner_index_c,values_c,total_values_c)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices computed for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_b
    !> This is the total number of values in C.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    !> This is the total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

    INCLUDE "sparse_includes/PairwiseMultiplyVectors.f90"

  END SUBROUTINE PairwiseMultiplyVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SVectorModule
