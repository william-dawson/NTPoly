!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the sparse matrix data type.
MODULE SparseMatrixModule_wrp
  USE SparseMatrixModule, ONLY : SparseMatrix_t, &
       & ConstructEmptySparseMatrix, ConstructSparseMatrixFromFile, &
       & ConstructFromTripletList, DestructSparseMatrix, CopySparseMatrix, &
       & ScaleSparseMatrix, IncrementSparseMatrix, Gemm, &
       & TransposeSparseMatrix, PrintSparseMatrix, MatrixToTripletList, &
       & GetRows, GetColumns, ComposeSparseMatrixColumns, DotSparseMatrix, &
       & PairwiseMultiplySparseMatrix
  USE TripletListModule_wrp, ONLY : TripletList_wrp
  USE MatrixMemoryPoolModule_wrp, ONLY : MatrixMemoryPool_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE WrapperModule, ONLY : SIZE_wrp
  USE iso_c_binding, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: SparseMatrix_wrp
     !> Actual data.
     TYPE(SparseMatrix_t), POINTER :: data
  END TYPE SparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptySparseMatrix_wrp
  PUBLIC :: ConstructSparseMatrixFromFile_wrp
  PUBLIC :: ConstructFromTripletList_wrp
  PUBLIC :: DestructSparseMatrix_wrp
  PUBLIC :: GetRows_wrp
  PUBLIC :: GetColumns_wrp
  PUBLIC :: CopySparseMatrix_wrp
  PUBLIC :: ScaleSparseMatrix_wrp
  PUBLIC :: IncrementSparseMatrix_wrp
  PUBLIC :: PairwiseMultiplySparseMatrix_wrp
  PUBLIC :: DotSparseMatrix_wrp
  PUBLIC :: Gemm_wrp
  PUBLIC :: TransposeSparseMatrix_wrp
  PUBLIC :: PrintSparseMatrix_wrp
  PUBLIC :: PrintSparseMatrixF_wrp
  PUBLIC :: MatrixToTripletList_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the empty sparse matrix constructor.
  !! @param[out] ih_this handle to the matrix being created.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmptySparseMatrix_wrp(ih_this, columns, rows) &
       & bind(c,name="ConstructEmptySparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: columns
    INTEGER(kind=c_int), INTENT(in) :: rows
    TYPE(SparseMatrix_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructEmptySparseMatrix(h_this%data,columns,rows)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructEmptySparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] ih_this the matrix being constructed.
  !! @param[in] file_name name of the file.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE ConstructSparseMatrixFromFile_wrp(ih_this, file_name, name_size) &
       & bind(c,name="ConstructSparseMatrixFromFile_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(SparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructSparseMatrixFromFile(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructSparseMatrixFromFile_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! @param[inout] ih_this handle to the matrix being constructed
  !! @param[in] ih_triplet_list handle to a sorted list of triplet values
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructFromTripletList_wrp(ih_this, ih_triplet_list, &
       & rows, columns) bind(c,name="ConstructFromTripletList_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_triplet_list(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: columns
    INTEGER(kind=c_int), INTENT(in) :: rows
    !! Local Data
    TYPE(SparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp)  :: h_triplet_list

    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    ALLOCATE(h_this%data)
    CALL ConstructFromTripletList(h_this%data, h_triplet_list%data, &
         & rows, columns)
    h_this = TRANSFER(ih_this,h_this)
  END SUBROUTINE ConstructFromTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix
  !! @param[inout] ih_this handle to the matrix to free up.
  PURE SUBROUTINE DestructSparseMatrix_wrp(ih_this) &
       & bind(c,name="DestructSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructSparseMatrix(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the copy a sparse matrix routine.
  !! @param[in] ih_matA handle to the matrix to copy
  !! @param[inout] ih_matB matB = matA
  PURE SUBROUTINE CopySparseMatrix_wrp(ih_matA, ih_matB) &
       & bind(c,name="CopySparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopySparseMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopySparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the row accessor.
  !! @param[in] ih_this handle to the matrix.
  !! @param[out] rows the number of rows.
  PURE SUBROUTINE GetRows_wrp(ih_this, rows) &
       & bind(c,name="GetRows_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out) :: rows
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    rows = GetRows(h_this%data)
  END SUBROUTINE GetRows_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the column accessor.
  !! @param[in] ih_this handle to the matrix.
  !! @param[out] columns the number of columns.
  PURE SUBROUTINE GetColumns_wrp(ih_this, columns) &
       & bind(c,name="GetColumns_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out) :: columns
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    columns = GetColumns(h_this%data)
  END SUBROUTINE GetColumns_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  !! @param[inout] ih_this handle to the matrix to scale.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleSparseMatrix_wrp(ih_this, constant) &
       & bind(c,name="ScaleSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: constant
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleSparseMatrix(h_this%data,constant)
  END SUBROUTINE ScaleSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in,out] ih_matB handle to Matrix B.
  !! @param[in] alpha_in multiplier.
  !! @param[in] threshold_in for flushing values to zero.
  PURE SUBROUTINE IncrementSparseMatrix_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) bind(c,name="IncrementSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: alpha_in
    REAL(NTREAL), INTENT(in) :: threshold_in
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementSparseMatrix(h_matA%data, h_matB%data, alpha_in, threshold_in)
  END SUBROUTINE IncrementSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in,out] ih_matB handle to Matrix B.
  PURE FUNCTION DotSparseMatrix_wrp(ih_matA, ih_matB) RESULT(product) &
       & bind(c,name="DotSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotSparseMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in] ih_matB handle to Matrix B.
  !! @param[out] ih_matC = A pairwise B
  SUBROUTINE PairwiseMultiplySparseMatrix_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="PairwiseMultiplySparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB
    TYPE(SparseMatrix_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplySparseMatrix(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE PairwiseMultiplySparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in] ih_matB handle to Matrix B.
  !! @param[out] ih_matC = alpha*matA*op( matB ) + beta*matC.
  !! @param[in] IsATransposed true if A is already transposed.
  !! @param[in] IsBTransposed true if B is already transposed.
  !! @param[in] alpha scales the multiplication.
  !! @param[in] beta scales matrix we sum on to.
  !! @param[in] threshold for flushing values to zero.
  !! @param[inout] ih_blocked_memory_pool handle to memory pool.
  SUBROUTINE Gemm_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & bind(c,name="Gemm_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(in) :: IsATransposed
    LOGICAL(kind=c_bool), INTENT(in) :: IsBTransposed
    REAL(NTREAL), INTENT(in) :: alpha
    REAL(NTREAL), INTENT(in) :: beta
    REAL(NTREAL), INTENT(in) :: threshold
    INTEGER(kind=c_int), INTENT(inout) :: ih_blocked_memory_pool(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB
    TYPE(SparseMatrix_wrp) :: h_matC
    TYPE(MatrixMemoryPool_wrp) :: h_blocked_memory_pool

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_blocked_memory_pool = TRANSFER(ih_blocked_memory_pool, &
         & h_blocked_memory_pool)

    CALL Gemm(h_matA%data, h_matB%data, h_matC%data, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%data)
  END SUBROUTINE Gemm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the matrix transpose function.
  !! @param[in] ih_matA handle to the matrix to be transposed.
  !! @param[inout] ih_matAT handle to the the input matrix transposed.
  PURE SUBROUTINE TransposeSparseMatrix_wrp(ih_matA, ih_matAT) &
       & bind(c,name="TransposeSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matAT(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matAT

    h_matA  = TRANSFER(ih_matA,h_matA)
    h_matAT = TRANSFER(ih_matAT,h_matAT)
    CALL TransposeSparseMatrix(h_matA%data,h_matAT%data)
  END SUBROUTINE TransposeSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints out a sparse matrix to file.
  !! @param[in] ih_this the matrix to be printed.
  !! @param[in] file_name optionally, you can pass a file to print to.
  !! @param[in] name_size the number of characters in the file name.
  SUBROUTINE PrintSparseMatrixF_wrp(ih_this, file_name, name_size) &
       & bind(c,name="PrintSparseMatrixF_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(in) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(in) :: name_size
    TYPE(SparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintSparseMatrix(h_this%data,local_string)
  END SUBROUTINE PrintSparseMatrixF_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints the sparse matrix to the console.
  !! @param[in] ih_this the matrix to be printed.
  SUBROUTINE PrintSparseMatrix_wrp(ih_this) &
       & bind(c,name="PrintSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintSparseMatrix(h_this%data)
  END SUBROUTINE PrintSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the routine that constructs a triplet list from a matrix.
  !! @param[in] ih_this handle to the matrix to construct the triplet list from
  !! @param[out] ih_triplet_list handle to the triplet list we created.
  SUBROUTINE MatrixToTripletList_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="MatrixToTripletList_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(out)   :: ih_triplet_list(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp)  :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_triplet_list%data)

    CALL MatrixToTripletList(h_this%data,h_triplet_list%data)

    ih_triplet_list = TRANSFER(ih_triplet_list,ih_triplet_list)
  END SUBROUTINE MatrixToTripletList_wrp
END MODULE SparseMatrixModule_wrp
