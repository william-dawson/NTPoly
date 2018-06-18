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

  END FUNCTION ConstructMatrixMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  PURE SUBROUTINE DestructMatrixMemoryPool(this)
    !! Parameters
    TYPE(MPOOLTYPE), INTENT(INOUT) :: this


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



  END FUNCTION CheckMemoryPoolValidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetPoolSparsity(this,sparsity)
    !! Parameters
    TYPE(MPOOLTYPE), INTENT(INOUT), TARGET :: this
    REAL(NTREAL), INTENT(IN) :: sparsity
    
  END SUBROUTINE SetPoolSparsity
