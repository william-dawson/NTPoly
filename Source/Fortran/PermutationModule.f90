!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for load balancing the matrix multiplication calculation.
MODULE PermutationModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure for storing permutations.
  TYPE, PUBLIC :: Permutation_t
     !> For each row/column, what index does it correspond to in the
     !! unperturbed matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: index_lookup
     !> For each row/column in the unperturbed, what index does it correspond to
     !! in this matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: reverse_index_lookup
  END TYPE Permutation_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDefaultPermutation
  PUBLIC :: ConstructReversePermutation
  PUBLIC :: ConstructRandomPermutation
  PUBLIC :: ConstructLimitedRandomPermutation
  PUBLIC :: DestructPermutation
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that preserves the original order.
  !! @param[inout] this the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  SUBROUTINE ConstructDefaultPermutation(this, matrix_dimension)
    !! Parameters
    TYPE(Permutation_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: matrix_dimension
    !! Local Data
    INTEGER :: counter

    ALLOCATE(this%index_lookup(matrix_dimension))
    ALLOCATE(this%reverse_index_lookup(matrix_dimension))

    !! Fill by counting.
    fill: DO counter = 1, matrix_dimension
       this%index_lookup(counter) = counter
       this%reverse_index_lookup(counter) = counter
    END DO fill

  END SUBROUTINE ConstructDefaultPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that reverses the original order.
  !! @param[inout] this the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  SUBROUTINE ConstructReversePermutation(this, matrix_dimension)
    !! Parameters
    TYPE(Permutation_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: matrix_dimension
    !! Local Data
    INTEGER :: counter

    ALLOCATE(this%index_lookup(matrix_dimension))
    ALLOCATE(this%reverse_index_lookup(matrix_dimension))

    !! Fill by counting.
    fill: DO counter = 1, matrix_dimension
       this%index_lookup(counter) = matrix_dimension - counter + 1
       this%reverse_index_lookup(counter) = counter
    END DO fill

  END SUBROUTINE ConstructReversePermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that has a random order.
  !! @param[inout] this the permutation to construct.
  !! @param[in] matrix_dimension size of the matrix.
  !! Implements Knuth shuffle.
  SUBROUTINE ConstructRandomPermutation(this, matrix_dimension)
    !! Parameters
    TYPE(Permutation_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: matrix_dimension
    !! Local Data
    INTEGER :: counter
    INTEGER :: random_integer
    REAL(KIND=8) :: rand_temp
    INTEGER :: swap_space
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    INTEGER :: seed_size

    !! Temporary, seed the random number generator
    CALL RANDOM_SEED(size=seed_size)
    ALLOCATE(seed(seed_size))
    seed = 0
    CALL RANDOM_SEED(put=seed)

    !! First fill by counting.
    CALL ConstructDefaultPermutation(this,matrix_dimension)

    !! Do the shuffle
    shuffle: DO counter=matrix_dimension,1,-1
       CALL RANDOM_NUMBER(rand_temp)
       random_integer = FLOOR(matrix_dimension*rand_temp)+1
       swap_space = this%index_lookup(matrix_dimension)
       this%index_lookup(matrix_dimension) = this%index_lookup(random_integer)
       this%index_lookup(random_integer) = swap_space
    END DO shuffle

    !! Compute the reverse lookup
    reverse: DO counter=1,matrix_dimension
       this%reverse_index_lookup(this%index_lookup(counter)) = counter
    END DO reverse
  END SUBROUTINE ConstructRandomPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that has a random order, but there is no
  !! permutation from beyond the actual matrix dimension.
  !! @param[inout] this the permutation to construct.
  !! @param[in] actual_matrix_dimension actual size of the matrix.
  !! @param[in] logical_matrix_dimension padded size of the matrix.
  SUBROUTINE ConstructLimitedRandomPermutation(this, actual_matrix_dimension, &
       & logical_matrix_dimension)
    !! Parameters
    TYPE(Permutation_t), INTENT(inout) :: this
    INTEGER, INTENT(in) :: actual_matrix_dimension
    INTEGER, INTENT(in) :: logical_matrix_dimension
    !! Local Data
    INTEGER :: counter
    INTEGER :: random_integer
    REAL(KIND=8) :: rand_temp
    INTEGER :: swap_space
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    INTEGER :: seed_size

    !! Temporary, seed the random number generator
    CALL RANDOM_SEED(size=seed_size)
    ALLOCATE(seed(seed_size))
    seed = 0
    CALL RANDOM_SEED(put=seed)

    !! First fill by counting.
    CALL ConstructDefaultPermutation(this,logical_matrix_dimension)

    !! Do the shuffle
    shuffle: DO counter=actual_matrix_dimension,1,-1
       CALL RANDOM_NUMBER(rand_temp)
       random_integer = FLOOR(actual_matrix_dimension*rand_temp)+1
       swap_space = this%index_lookup(actual_matrix_dimension)
       this%index_lookup(actual_matrix_dimension) = &
            & this%index_lookup(random_integer)
       this%index_lookup(random_integer) = swap_space
    END DO shuffle

    !! Compute the reverse lookup
    reverse: DO counter=1,logical_matrix_dimension
       this%reverse_index_lookup(this%index_lookup(counter)) = counter
    END DO reverse
  END SUBROUTINE ConstructLimitedRandomPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  !! @param[inout] this the permutation to destruct.
  PURE SUBROUTINE DestructPermutation(this)
    !! Parameters
    TYPE(Permutation_t), INTENT(inout) :: this

    IF (ALLOCATED(this%index_lookup)) THEN
       DEALLOCATE(this%index_lookup)
    END IF
    IF (ALLOCATED(this%reverse_index_lookup)) THEN
       DEALLOCATE(this%reverse_index_lookup)
    END IF
  END SUBROUTINE DestructPermutation
END MODULE PermutationModule
