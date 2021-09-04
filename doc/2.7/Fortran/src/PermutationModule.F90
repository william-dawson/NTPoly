!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for load balancing the matrix multiplication calculation.
MODULE PermutationModule
  USE DataTypesModule, ONLY : NTREAL
  USE ProcessGridModule, ONLY : global_grid, ProcessGrid_t
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure for storing permutations.
  TYPE, PUBLIC :: Permutation_t
     !> For each row/column, what index does it correspond to in the
     !> unperturbed matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: index_lookup
     !> For each row/column in the unperturbed, what index does it correspond to
     !> in this matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: reverse_index_lookup
  END TYPE Permutation_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDefaultPermutation
  PUBLIC :: ConstructReversePermutation
  PUBLIC :: ConstructRandomPermutation
  PUBLIC :: ConstructLimitedRandomPermutation
  PUBLIC :: DestructPermutation
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that preserves the original order.
  SUBROUTINE ConstructDefaultPermutation(this, matrix_dimension)
    !> The permutation to construct.
    TYPE(Permutation_t), INTENT(INOUT) :: this
    !> The dimension of the matrix.
    INTEGER, INTENT(IN) :: matrix_dimension
    !! Local Data
    INTEGER :: counter

    CALL DestructPermutation(this)

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
  SUBROUTINE ConstructReversePermutation(this, matrix_dimension)
    !> A permutation that reverses the original order.
    TYPE(Permutation_t), INTENT(INOUT) :: this
    !> The size of the matrix.
    INTEGER, INTENT(IN) :: matrix_dimension
    !! Local Data
    INTEGER :: counter

    CALL DestructPermutation(this)

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
  !> Implements Knuth shuffle.
  SUBROUTINE ConstructRandomPermutation(this, matrix_dimension, &
       & process_grid_in)
    !> A permutation that reverses the original order.
    TYPE(Permutation_t), INTENT(INOUT) :: this
    !> The size of the matrix.
    INTEGER, INTENT(IN) :: matrix_dimension
    !> A permutation should be shared amongst these processes.
    !> This is to synchronize random number across processes.
    TYPE(ProcessGrid_t), INTENT(INOUT), OPTIONAL :: process_grid_in
    !! Local Data
    INTEGER :: counter
    INTEGER :: random_integer
    REAL(KIND=NTREAL) :: rand_temp
    INTEGER :: swap_space
    INTEGER :: ierr

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

    !! Broadcast the lookup (so each process has the same value)
    IF (PRESENT(process_grid_in)) THEN
       CALL MPI_Bcast(this%index_lookup, SIZE(this%index_lookup), &
            & MPI_INTEGER, 0, process_grid_in%global_comm, ierr)
    ELSE
       CALL MPI_Bcast(this%index_lookup, SIZE(this%index_lookup), &
            & MPI_INTEGER, 0, global_grid%global_comm, ierr)
    END IF

    !! Compute the reverse lookup
    reverse: DO counter=1,matrix_dimension
       this%reverse_index_lookup(this%index_lookup(counter)) = counter
    END DO reverse
  END SUBROUTINE ConstructRandomPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Constructs a permutation that has a random order, but there is no
  !> permutation from beyond the actual matrix dimension.
  SUBROUTINE ConstructLimitedRandomPermutation(this, actual_matrix_dimension, &
       & logical_matrix_dimension, process_grid_in)
    !> The permutation to construct.
    TYPE(Permutation_t), INTENT(inout) :: this
    !> Actual size of the matrix.
    INTEGER, INTENT(IN) :: actual_matrix_dimension
    !> Padded size of the matrix.
    INTEGER, INTENT(IN) :: logical_matrix_dimension
    !> A permutation should be shared amongst these processes.
    !> This is to synchronize random number across processes.
    TYPE(ProcessGrid_t), INTENT(INOUT), OPTIONAL :: process_grid_in
    !! Local Data
    INTEGER :: counter
    INTEGER :: random_integer
    REAL(KIND=NTREAL) :: rand_temp
    INTEGER :: swap_space
    INTEGER :: ierr

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

    !! Broadcast the lookup (so each process has the same value)
    IF (PRESENT(process_grid_in)) THEN
       CALL MPI_Bcast(this%index_lookup, SIZE(this%index_lookup), &
            & MPI_INTEGER, 0, process_grid_in%global_comm, ierr)
    ELSE
       CALL MPI_Bcast(this%index_lookup, SIZE(this%index_lookup), &
            & MPI_INTEGER, 0, global_grid%global_comm, ierr)
    END IF

    !! Compute the reverse lookup
    reverse: DO counter=1,logical_matrix_dimension
       this%reverse_index_lookup(this%index_lookup(counter)) = counter
    END DO reverse
  END SUBROUTINE ConstructLimitedRandomPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  PURE SUBROUTINE DestructPermutation(this)
    !> The permutation to destruct.
    TYPE(Permutation_t), INTENT(inout) :: this

    IF (ALLOCATED(this%index_lookup)) THEN
       DEALLOCATE(this%index_lookup)
    END IF
    IF (ALLOCATED(this%reverse_index_lookup)) THEN
       DEALLOCATE(this%reverse_index_lookup)
    END IF
  END SUBROUTINE DestructPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PermutationModule
