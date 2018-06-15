!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @param[in] sparsity_in estimated sparsity (optional).
  FUNCTION ConstructMatrixMemoryPool(columns, rows, sparsity_in) RESULT(this)
    !! Parameters
    TYPE(MPOOLTYPE), TARGET :: this
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in
    !! Temporary variables
    INTEGER :: alloc_stat
    INTEGER :: num_buckets

    this%columns = columns
    this%rows = rows

    IF (.NOT. PRESENT(sparsity_in)) THEN
       this%hash_size = 1
    ELSE
       this%hash_size = INT(1.0/sparsity_in)
       IF (this%hash_size > columns) this%hash_size = columns
    END IF

    num_buckets = columns/this%hash_size + 1

    !! Allocate
    ALLOCATE(this%pruned_list(columns*rows), stat=alloc_stat)
    ALLOCATE(this%value_array(columns,rows), stat=alloc_stat)
    ALLOCATE(this%dirty_array(columns,rows), stat=alloc_stat)

    ALLOCATE(this%hash_index(columns,rows))
    ALLOCATE(this%inserted_per_bucket(columns,rows))

    this%value_array = 0
    this%hash_index = 0
    this%inserted_per_bucket = 0
    this%dirty_array = .FALSE.
  END FUNCTION ConstructMatrixMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  PURE SUBROUTINE DestructMatrixMemoryPool(this)
    !! Parameters
    TYPE(MPOOLTYPE), INTENT(INOUT) :: this

    !! Perform deallocations.
    IF (ALLOCATED(this%pruned_list)) DEALLOCATE(this%pruned_list)
    IF (ALLOCATED(this%value_array)) DEALLOCATE(this%value_array)
    IF (ALLOCATED(this%dirty_array)) DEALLOCATE(this%dirty_array)
    IF (ALLOCATED(this%hash_index)) DEALLOCATE(this%hash_index)
    IF (ALLOCATED(this%inserted_per_bucket)) &
         & DEALLOCATE(this%inserted_per_bucket)
  END SUBROUTINE DestructMatrixMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckMemoryPoolValidity(this, columns, rows) RESULT(isvalid)
    !! Parameters
    TYPE(MPOOLTYPE), INTENT(in) :: this
    INTEGER, INTENT(IN) :: columns
    INTEGER, INTENT(IN) :: rows
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (.NOT. ALLOCATED(this%pruned_list)) isvalid = .FALSE.
    IF (.NOT. ALLOCATED(this%value_array)) isvalid = .FALSE.

    !! Check allocation size
    IF (.NOT. SIZE(this%value_array,dim=2) .EQ. rows) THEN
       isvalid = .FALSE.
    END IF
    IF (.NOT. SIZE(this%value_array,dim=1) .EQ. columns) THEN
       isvalid = .FALSE.
    END IF

  END FUNCTION CheckMemoryPoolValidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetPoolSparsity(this,sparsity)
    !! Parameters
    TYPE(MPOOLTYPE), INTENT(INOUT), TARGET :: this
    REAL(NTREAL), INTENT(IN) :: sparsity
    !! Local Variables
    INTEGER :: num_buckets

    this%hash_size = INT(1.0/sparsity)
    IF (this%hash_size > this%columns) this%hash_size = this%columns
    num_buckets = this%columns/this%hash_size + 1
  END SUBROUTINE SetPoolSparsity
