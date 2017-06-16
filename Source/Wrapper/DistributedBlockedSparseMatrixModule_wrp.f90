!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Blocked Sparse Matrix.
MODULE DistributedBlockedSparseMatrixModule_wrp
#ifdef NOBLOCK
  USE DistributedSparseMatrixModule
#else
  USE DistributedBlockedSparseMatrixModule
#endif
  USE DistributedMatrixMemoryPoolModule_wrp, ONLY : &
       & DistributedMatrixMemoryPool_wrp
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE TripletListModule_wrp, ONLY : TripletList_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: DistributedSparseMatrix_wrp
     !> Actual data.
     TYPE(DistributedSparseMatrix), POINTER :: data
  END TYPE DistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmpty_wrp
  PUBLIC :: CopyDistributedSparseMatrix_wrp
  PUBLIC :: DestructDistributedSparseMatrix_wrp
  PUBLIC :: ConstructFromMatrixMarket_wrp
  PUBLIC :: ConstructFromBinary_wrp
  PUBLIC :: WriteToBinary_wrp
  PUBLIC :: WriteToMatrixMarket_wrp
  PUBLIC :: FillFromTripletList_wrp
  PUBLIC :: FillDistributedPermutation_wrp
  PUBLIC :: FillDistributedIdentity_wrp
  PUBLIC :: GetActualDimension_wrp
  PUBLIC :: GetLogicalDimension_wrp
  PUBLIC :: GetTripletList_wrp
  PUBLIC :: IncrementDistributedSparseMatrix_wrp
  PUBLIC :: DotDistributedSparseMatrix_wrp
  PUBLIC :: DistributedPairwiseMultiply_wrp
  PUBLIC :: DistributedGemm_wrp
  PUBLIC :: ScaleDistributedSparseMatrix_wrp
  PUBLIC :: DistributedSparseNorm_wrp
  PUBLIC :: Trace_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the constructor of an empty sparse, distributed, matrix.
  !! @param[inout] ih_this the matrix to be constructed.
  !! @param[in] matrix_dim the dimension of the full matrix.
  PURE SUBROUTINE ConstructEmpty_wrp(ih_this,matrix_dim) &
       & bind(c,name="ConstructEmpty_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: matrix_dim
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructEmpty(h_this%data,matrix_dim)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructEmpty_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  !! @param[out] ih_this the file being constructed.
  !! @param[in] file_name name of the file to read.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE ConstructFromMatrixMarket_wrp(ih_this,file_name,name_size) &
       & bind(c,name="ConstructFromMatrixMarket_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructFromMatrixMarket(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructFromMatrixMarket_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  !! @param[out] ih_this the file being constructed.
  !! @param[in] file_name name of the file to read.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE ConstructFromBinary_wrp(ih_this,file_name,name_size) &
       & bind(c,name="ConstructFromBinary_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructFromBinary(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructFromBinary_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  !! @param[in] ih_matA matrix to copy
  !! @param[inout] ih_matB = matA
  PURE SUBROUTINE CopyDistributedSparseMatrix_wrp(ih_matA,ih_matB) &
       & bind(c,name="CopyDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopyDistributedSparseMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopyDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix
  !! @param[in,out] ih_this the matrix to destruct
  PURE SUBROUTINE DestructDistributedSparseMatrix_wrp(ih_this) &
       & bind(c,name="DestructDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructDistributedSparseMatrix(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a file.
  !! @param[in] ih_this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE WriteToBinary_wrp(ih_this,file_name,name_size) &
       & bind(c,name="WriteToBinary_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteToBinary(h_this%data,local_string)
  END SUBROUTINE WriteToBinary_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a matrix market file.
  !! @param[in] ih_this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE WriteToMatrixMarket_wrp(ih_this,file_name,name_size) &
       & bind(c,name="WriteToMatrixMarket_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteToMatrixMarket(h_this%data,local_string)
  END SUBROUTINE WriteToMatrixMarket_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists.
  !! @param[inout] ih_this the matrix being filled.
  !! @param[in] ih_triplet_list the triplet list of values.
  SUBROUTINE FillFromTripletList_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="FillFromTripletList_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_triplet_list(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL FillFromTripletList(h_this%data, h_triplet_list%data)
  END SUBROUTINE FillFromTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] ih_this the matrix being filled.
  PURE SUBROUTINE FillDistributedIdentity_wrp(ih_this) &
       & bind(c,name="FillDistributedIdentity_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL FillDistributedIdentity(h_this%data)
  END SUBROUTINE FillDistributedIdentity_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation matrix.
  !! @param[inout] ih_this the matrix being filled.
  !! @param[in] ih_permutation the perumutation to fill with.
  !! @param[in] permute_rows whether to permute rows or columns.
  SUBROUTINE FillDistributedPermutation_wrp(ih_this, ih_permutation, &
       & permute_rows) bind(c,name="FillDistributedPermutation_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_permutation(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(in) :: permute_rows
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_permutation

    h_this = TRANSFER(ih_this,h_this)
    h_permutation = TRANSFER(ih_permutation,h_permutation)

    CALL FillDistributedPermutation(h_this%data,h_permutation%data%index_lookup, &
         & permuterows = LOGICAL(permute_rows))
  END SUBROUTINE FillDistributedPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Actual Dimension accessor.
  !! @param[in] ih_this handle to the matrix.
  !! @param[out] mat_dimension the size of the matrix.
  PURE SUBROUTINE GetActualDimension_wrp(ih_this, mat_dimension) &
       & bind(c,name="GetActualDimension_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out) :: mat_dimension
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetActualDimension(h_this%data)
  END SUBROUTINE GetActualDimension_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Logical Dimension accessor.
  !! @param[in] ih_this handle to the matrix.
  !! @param[out] mat_dimension the size of the matrix.
  PURE SUBROUTINE GetLogicalDimension_wrp(ih_this, mat_dimension) &
       & bind(c,name="GetLogicalDimension_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out) :: mat_dimension
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetLogicalDimension(h_this%data)
  END SUBROUTINE GetLogicalDimension_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !! @param[in] ih_this the distributed sparse matrix.
  !! @param[inout] ih_triplet_list the list to fill.
  PURE SUBROUTINE GetTripletList_wrp(ih_this, ih_triplet_list)
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_triplet_list(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetTripletList(h_this%data,h_triplet_list%data)
  END SUBROUTINE GetTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !! @param[in] ih_matA Matrix A
  !! @param[in,out] ih_matB Matrix B
  !! @param[in] alpha_in multiplier.
  !! @param[in] threshold_in for flushing values to zero.
  SUBROUTINE IncrementDistributedSparseMatrix_wrp(ih_matA, ih_matB,&
       & alpha_in,threshold_in) &
       & bind(c,name="IncrementDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: alpha_in
    REAL(NTREAL), INTENT(in) :: threshold_in
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementDistributedSparseMatrix(h_matA%data, h_matB%data, alpha_in, &
         & threshold_in)
  END SUBROUTINE IncrementDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(matA,matB)
  !! @param[in] ih_matA Matrix A
  !! @param[in] ih_matB Matrix B
  !! @result product dot product of the two matrices
  FUNCTION DotDistributedSparseMatrix_wrp(ih_matA, ih_matB) RESULT(product)&
       & bind(c,name="DotDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotDistributedSparseMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise multiply two matrices together.
  !! @param[in] ih_matA Matrix A
  !! @param[in] ih_matB Matrix B
  !! @param[out] ih_matC = matA*matB
  SUBROUTINE DistributedPairwiseMultiply_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="DistributedPairwiseMultiply_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB
    TYPE(DistributedSparseMatrix_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL DistributedPairwiseMultiply(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE DistributedPairwiseMultiply_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! @param[in] ih_matA Matrix A
  !! @param[in] ih_matB Matrix B
  !! @param[out] ih_matC = alpha*matA*matB + beta*matC
  !! @param[in] alpha_in scales the multiplication
  !! @param[in] beta_in scales matrix we sum on to
  !! @param[in] threshold_in for flushing values to zero.
  !! @param[inout] ih_memory_pool_in a memory pool that can be used for the
  !! calculation.
  SUBROUTINE DistributedGemm_wrp(ih_matA, ih_matB, ih_matC, alpha_in, &
       & beta_in, threshold_in, ih_memory_pool_in) &
       & bind(c,name="DistributedGemm_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: alpha_in
    REAL(NTREAL), INTENT(in) :: beta_in
    REAL(NTREAL), INTENT(in) :: threshold_in
    INTEGER(kind=c_int), INTENT(inout) :: ih_memory_pool_in(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB
    TYPE(DistributedSparseMatrix_wrp) :: h_matC
    TYPE(DistributedMatrixMemoryPool_wrp) :: h_memory_pool_in

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_memory_pool_in = TRANSFER(ih_memory_pool_in,h_memory_pool_in)

    CALL DistributedGemm(h_matA%data, h_matB%data, h_matC%data, &
         & alpha_in, beta_in, threshold_in, h_memory_pool_in%data)
  END SUBROUTINE DistributedGemm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  !! @param[inout] ih_this Matrix to scale.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleDistributedSparseMatrix_wrp(ih_this, constant) &
       & bind(c,name="ScaleDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: constant
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleDistributedSparseMatrix(h_this%data,constant)
  END SUBROUTINE ScaleDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  !! @param[in] ih_this the matrix to compute the norm of.
  !! @return the norm value of the full distributed sparse matrix.
  FUNCTION DistributedSparseNorm_wrp(ih_this) &
       & bind(c,name="DistributedSparseNorm_wrp") RESULT(norm_value)
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: norm_value
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    norm_value = DistributedSparseNorm(h_this%data)
  END FUNCTION DistributedSparseNorm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  !! @param[in] ih_this the matrix to compute the norm of.
  !! @return the trace value of the full distributed sparse matrix.
  FUNCTION Trace_wrp(ih_this) &
       & bind(c,name="Trace_wrp") RESULT(trace_value)
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: trace_value
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    trace_value = Trace(h_this%data)
  END FUNCTION Trace_wrp
END MODULE DistributedBlockedSparseMatrixModule_wrp
