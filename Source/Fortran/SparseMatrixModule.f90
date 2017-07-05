!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE SparseMatrixModule
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixMarketModule
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_t, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE SparseVectorModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY: TripletList_t, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletModule, ONLY : Triplet_t
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a CSR matrix.
  TYPE, PUBLIC :: SparseMatrix_t
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE SparseMatrix_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Construct/Destruct
  PUBLIC :: ConstructEmptySparseMatrix
  PUBLIC :: ConstructSparseMatrixFromFile
  PUBLIC :: ConstructFromTripletList
  PUBLIC :: DestructSparseMatrix
  PUBLIC :: CopySparseMatrix
  !! Basic Accessors
  PUBLIC :: GetRows
  PUBLIC :: GetColumns
  !! Linear Algebra
  PUBLIC :: ScaleSparseMatrix
  PUBLIC :: IncrementSparseMatrix
  PUBLIC :: DotSparseMatrix
  PUBLIC :: PairwiseMultiplySparseMatrix
  PUBLIC :: Gemm
  !public :: Gemm2
  PUBLIC :: SparseMatrixNorm
  PUBLIC :: SparseMatrixGrandSum
  !! Helper routines
  PUBLIC :: SplitSparseMatrixColumns
  PUBLIC :: ComposeSparseMatrixColumns
  PUBLIC :: TransposeSparseMatrix
  PUBLIC :: PrintSparseMatrix
  PUBLIC :: MatrixToTripletList
  PUBLIC :: CheckIfIdentity
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix with a certain number of columns and rows.
  !! Will allocate storage for the outer values.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmptySparseMatrix(this, columns, rows)
    !! Parameters
    TYPE(SparseMatrix_t),INTENT(out)  :: this
    INTEGER, INTENT(in) :: columns, rows

    CALL DestructSparseMatrix(this)
    this%rows = rows
    this%columns = columns
    ALLOCATE(this%outer_index(this%columns+1))
    this%outer_index = 0
  END SUBROUTINE ConstructEmptySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  SUBROUTINE ConstructSparseMatrixFromFile(this, file_name)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(out) :: this
    CHARACTER(len=*), INTENT(in)   :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    INTEGER :: temp_rows, temp_columns, temp_total_values
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    CHARACTER(len=81) :: input_buffer
    INTEGER :: file_handler
    TYPE(Triplet_t) :: temporary
    INTEGER :: counter
    LOGICAL :: found_comment_line
    LOGICAL :: error_occured
    file_handler = 16

    OPEN(file_handler,file=file_name,status='old')

    !! Parse the header.
    READ(file_handler,fmt='(A)') input_buffer
    error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
         & pattern_type)

    !! Extra Comment Lines
    found_comment_line = .TRUE.
    DO WHILE(found_comment_line)
       !read(file_handler,*) input_buffer
       READ(file_handler,fmt='(A)') input_buffer
       IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
          found_comment_line = .FALSE.
       END IF
    END DO

    !! Main data
    READ(input_buffer,*) temp_rows, temp_columns, temp_total_values
    CALL ConstructTripletList(triplet_list)

    !! Read Values
    DO counter = 1, temp_total_values
       READ(file_handler,*) temporary%index_row, temporary%index_column, &
            & temporary%point_value
       CALL AppendToTripletList(triplet_list,temporary)
    END DO
    CLOSE(file_handler)

    CALL SymmetrizeTripletList(triplet_list, pattern_type)
    CALL SortTripletList(triplet_list,temp_columns,sorted_triplet_list)
    CALL ConstructFromTripletList(this, sorted_triplet_list, temp_rows, &
         & temp_columns)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE ConstructSparseMatrixFromFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructFromTripletList(this,triplet_list,rows,columns)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(out) :: this
    TYPE(TripletList_t), INTENT(in) :: triplet_list
    INTEGER, INTENT(in)             :: rows, columns
    INTEGER :: outer_array_ptr
    INTEGER :: values_counter

    IF (ALLOCATED(this%outer_index)) THEN
       CALL DestructSparseMatrix(this)
    END IF

    this%rows = rows
    this%columns = columns

    !! Allocate
    ALLOCATE(this%outer_index(this%columns+1))
    this%outer_index = 0
    ALLOCATE(this%inner_index(triplet_list%CurrentSize))
    ALLOCATE(this%values(triplet_list%CurrentSize))

    !! Insert Values
    outer_array_ptr = 1
    DO values_counter = 1, triplet_list%CurrentSize
       !! Moving on to the next column?
       DO WHILE(.NOT. triplet_list%data(values_counter)%index_column .EQ. &
            & outer_array_ptr)
          outer_array_ptr = outer_array_ptr + 1
          this%outer_index(outer_array_ptr+1) = this%outer_index(outer_array_ptr)
       END DO
       this%outer_index(outer_array_ptr+1)=this%outer_index(outer_array_ptr+1)+1
       !! Insert inner index and value
       this%inner_index(values_counter) = &
            & triplet_list%data(values_counter)%index_row
       this%values(values_counter) = &
            & triplet_list%data(values_counter)%point_value
    END DO

    !! Fill In The Rest Of The Outer Values
    DO outer_array_ptr = outer_array_ptr+2, this%columns+1
       this%outer_index(outer_array_ptr) = this%outer_index(outer_array_ptr-1)
    END DO

  END SUBROUTINE ConstructFromTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! This will always check if arrays are actually allocated, so you can feel
  !! free to destruct a matrix even if it has no data.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE DestructSparseMatrix(this)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(inout) :: this

    IF (ALLOCATED(this%outer_index)) THEN
       DEALLOCATE(this%outer_index)
    END IF
    IF (ALLOCATED(this%inner_index)) THEN
       DEALLOCATE(this%inner_index)
    END IF
    IF (ALLOCATED(this%values)) THEN
       DEALLOCATE(this%values)
    END IF
  END SUBROUTINE DestructSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopySparseMatrix(matA, matB)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: matA
    TYPE(SparseMatrix_t), INTENT(inout) :: matB

    CALL DestructSparseMatrix(matB)
    matB = matA
  END SUBROUTINE CopySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetRows(this) RESULT(rows)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    INTEGER :: rows
    rows = this%rows
  END FUNCTION GetRows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetColumns(this) RESULT(columns)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    INTEGER :: columns
    columns = this%columns
  END FUNCTION GetColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleSparseMatrix(matA,constant)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(inout) :: matA
    REAL(NTREAL), INTENT(in) :: constant

    matA%values = constant * matA%values
  END SUBROUTINE ScaleSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier. Optional, default is 1.0
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @todo I don't like this hack where I have to check if MatrixB is allocated.
  PURE SUBROUTINE IncrementSparseMatrix(matA, matB, alpha_in, threshold_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in)  :: matA
    TYPE(SparseMatrix_t), INTENT(inout) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    TYPE(SparseMatrix_t) :: matC
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

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

    size_of_a = matA%outer_index(matA%columns+1)

    !! Allocate sufficient space for matC
    CALL ConstructEmptySparseMatrix(matC, matA%columns, matA%rows)
    IF (ALLOCATED(matB%values)) THEN
       size_of_b = matB%outer_index(matB%columns+1)
       ALLOCATE(matC%inner_index(size_of_a+size_of_b))
       ALLOCATE(matC%values(size_of_a+size_of_b))
    ELSE
       ALLOCATE(matC%inner_index(size_of_a))
       ALLOCATE(matC%values(size_of_a))
    END IF

    !! Perform loops
    total_counter_a = 1
    total_counter_b = 1
    total_counter_c = 1
    DO outer_counter = 1, matA%columns
       !! Inner counters
       elements_per_inner_a = matA%outer_index(outer_counter+1) - &
            & matA%outer_index(outer_counter)
       elements_per_inner_b = matB%outer_index(outer_counter+1) - &
            & matB%outer_index(outer_counter)
       CALL AddSparseVectors(&
            matA%inner_index(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matA%values(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matB%inner_index(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matB%values(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matC%inner_index(total_counter_c:),matC%values(total_counter_c:),&
            indices_added_into_c, alpha, threshold)
       matC%outer_index(outer_counter+1) = matC%outer_index(outer_counter)+&
            & indices_added_into_c
       total_counter_a = total_counter_a + elements_per_inner_a
       total_counter_b = total_counter_b + elements_per_inner_b
       total_counter_c = total_counter_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructSparseMatrix(matB)
    CALL ConstructEmptySparseMatrix(matB,matC%columns,matC%rows)
    matB%outer_index = matC%outer_index
    ALLOCATE(matB%inner_index(matC%outer_index(matC%columns+1)))
    ALLOCATE(matB%values(matC%outer_index(matC%columns+1)))
    matB%inner_index = matC%inner_index(:matC%outer_index(matC%columns+1))
    matB%values = matC%values(:matC%outer_index(matC%columns+1))
    CALL DestructSparseMatrix(matC)
  END SUBROUTINE IncrementSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiplySparseMatrix(matA, matB, matC)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in)  :: matA
    TYPE(SparseMatrix_t), INTENT(in) :: matB
    TYPE(SparseMatrix_t), INTENT(inout) :: matC
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    TYPE(SparseMatrix_t) :: TempMat
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    CALL ConstructEmptySparseMatrix(TempMat, matA%columns, matA%rows)
    size_of_a = matA%outer_index(matA%columns+1)
    size_of_b = matB%outer_index(matB%columns+1)
    ALLOCATE(TempMat%inner_index(MIN(size_of_a,size_of_b)))
    ALLOCATE(TempMat%values(MIN(size_of_a,size_of_b)))

    !! Perform loops
    total_counter_a = 1
    total_counter_b = 1
    total_counter_c = 1
    DO outer_counter = 1, matA%columns
       !! Inner counters
       elements_per_inner_a = matA%outer_index(outer_counter+1) - &
            & matA%outer_index(outer_counter)
       elements_per_inner_b = matB%outer_index(outer_counter+1) - &
            & matB%outer_index(outer_counter)
       CALL PairwiseMultiplyVectors(&
            matA%inner_index(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matA%values(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matB%inner_index(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matB%values(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            TempMat%inner_index(total_counter_c:),TempMat%values(total_counter_c:),&
            indices_added_into_c)
       TempMat%outer_index(outer_counter+1) = TempMat%outer_index(outer_counter)+&
            & indices_added_into_c
       total_counter_a = total_counter_a + elements_per_inner_a
       total_counter_b = total_counter_b + elements_per_inner_b
       total_counter_c = total_counter_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructSparseMatrix(matC)
    CALL ConstructEmptySparseMatrix(matC,TempMat%columns,TempMat%rows)
    matC%outer_index = TempMat%outer_index
    ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
    ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
    matC%inner_index = TempMat%inner_index(:TempMat%outer_index(TempMat%columns+1))
    matC%values = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
    CALL DestructSparseMatrix(TempMat)
  END SUBROUTINE PairwiseMultiplySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @result product
  PURE FUNCTION DotSparseMatrix(matA, matB) RESULT(product)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: matA
    TYPE(SparseMatrix_t), INTENT(in) :: matB
    REAL(NTREAL) :: product
    !! Local Variables
    TYPE(SparseMatrix_t) :: matC

    CALL PairwiseMultiplySparseMatrix(matA,matB,matC)

    product = SparseMatrixGrandSum(matC)
    CALL DestructSparseMatrix(matC)
  END FUNCTION DotSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*op( matB ) + beta*matC
  !! This version is not really linear scaling, but it should be fast enough.
  !! Basically, we create a big buffer region of zeros so that we can accumulate
  !! in O(1) time. Then we scan over the buffer region and search for filled
  !! values. For small enough matrices this is fine, but it definitely isn't
  !! optimal.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[out] matC = alpha*matA*op( matB ) + beta*matC.
  !! @param[in] IsATransposed_in true if A is already transposed.
  !! @param[in] IsBTransposed_in true if B is already transposed.
  !! @param[in] alpha_in scales the multiplication.
  !! @param[in] beta_in scales matrix we sum on to.
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @param[inout] blocked_memory_pool_in an optional memory pool for doing the
  !! calculation.
  !! @todo more performance tuning.
  ! subroutine Gemm2(matA, matB, matC, IsATransposed_in, IsBTransposed_in, &
  !   & alpha_in, beta_in, threshold_in, blocked_memory_pool_in)
  !   !! Parameters
  !   type(SparseMatrix_t), intent(in)  :: matA
  !   type(SparseMatrix_t), intent(in)  :: matB
  !   type(SparseMatrix_t), intent(inout) :: matC
  !   logical, optional, intent(in) :: IsATransposed_in
  !   logical, optional, intent(in) :: IsBTransposed_in
  !   real(NTREAL), optional, intent(in) :: alpha_in
  !   real(NTREAL), optional, intent(in) :: beta_in
  !   real(NTREAL), optional, intent(in) :: threshold_in
  !   type(MatrixMemoryPool_t), optional, &
  !     & intent(inout), target :: blocked_memory_pool_in
  !   !! Intermediate Data
  !   type(SparseMatrix_t) :: matAB
  !   logical :: IsATransposed, IsBTransposed
  !   real(NTREAL) :: alpha
  !   real(NTREAL) :: beta
  !   real(NTREAL) :: threshold
  !   type(SparseMatrix_t) :: matAT, matBT
  !   type(MatrixMemoryPool_t) :: blocked_memory_pool
  !   !! Counters and temporary data
  !   integer :: mat_c_columns, mat_c_rows
  !   !! For Efficiency Purposes
  !   logical :: pool_flag
  !
  !   !! Process Optional Parameters
  !   if (.not. present(alpha_in)) then
  !     alpha = 1.0d+0
  !   else
  !     alpha = alpha_in
  !   end if
  !   if (.not. present(beta_in)) then
  !     beta = 0.0
  !   else
  !     beta = beta_in
  !   end if
  !   if (.not. present(IsATransposed_in)) then
  !     IsATransposed = .false.
  !   else
  !     IsATransposed = IsATransposed_in
  !   end if
  !   if (.not. present(IsBTransposed_in)) then
  !     IsBTransposed = .false.
  !   else
  !     IsBTransposed = IsBTransposed_in
  !   end if
  !   if (.not. present(threshold_in)) then
  !     threshold = 0.0
  !   else
  !     threshold = threshold_in
  !   end if
  !
  !   !! Storage details for result matrix
  !   if (IsATransposed) then
  !     mat_c_rows = matA%columns
  !   else
  !     mat_c_rows = matA%rows
  !   end if
  !   if (IsBTransposed) then
  !     mat_c_columns = matB%rows
  !   else
  !     mat_c_columns = matB%columns
  !   end if
  !
  !   !! Initialization of Memory
  !   !call StartTimer("Initialize Prune")
  !   if (.not. present(blocked_memory_pool_in)) then
  !     call ConstructMatrixMemoryPool(blocked_memory_pool,mat_c_columns, &
  !          & mat_c_rows)
  !     pool_flag = .false.
  !   elseif (.not. CheckMemoryPoolValidity(blocked_memory_pool_in, &
  !         & mat_c_columns, mat_c_rows)) then
  !     call DestructMatrixMemoryPool(blocked_memory_pool_in)
  !     call ConstructMatrixMemoryPool(blocked_memory_pool_in,mat_c_columns, &
  !          & mat_c_rows)
  !     pool_flag = .true.
  !   else
  !     pool_flag = .true.
  !   end if
  !
  !   !! Block A and B
  !   !call StartTimer("Transposes")
  !   if (.not. IsATransposed) then
  !     call TransposeSparseMatrix(matA,matAT)
  !   end if
  !   if (.not. IsBTransposed) then
  !     call TransposeSparseMatrix(matB,matBT)
  !   end if
  !   !call StopTimer("Transposes")
  !
  !   !call StopTimer("Initialize Prune")
  !
  !   !call StartTimer("MM Loop")
  !   !call StartTimer("Multiply Block")
  !   if (pool_flag) then
  !     if (IsATransposed .and. IsBTransposed) then
  !       call MultiplyBlock(matA, matB, blocked_memory_pool_in%value_array)
  !     elseif (IsATransposed) then
  !       call MultiplyBlock(matA, matBT, blocked_memory_pool_in%value_array)
  !     elseif (IsBTransposed) then
  !       call MultiplyBlock(matAT, matB, blocked_memory_pool_in%value_array)
  !     else
  !       call MultiplyBlock(matAT, matBT, blocked_memory_pool_in%value_array)
  !     end if
  !   else
  !     if (IsATransposed .and. IsBTransposed) then
  !       call MultiplyBlock(matA, matB, blocked_memory_pool%value_array)
  !     elseif (IsATransposed) then
  !       call MultiplyBlock(matA, matBT, blocked_memory_pool%value_array)
  !     elseif (IsBTransposed) then
  !       call MultiplyBlock(matAT, matB, blocked_memory_pool%value_array)
  !     else
  !       call MultiplyBlock(matAT, matBT, blocked_memory_pool%value_array)
  !     end if
  !   end if
  !   !call StopTimer("Multiply Block")
  !
  !   !! Go from triplets to return matrix
  !   !call StartTimer("Prune List")
  !   !call StartTimer("Prune List")
  !   if (pool_flag) then
  !     call PruneList(blocked_memory_pool_in,alpha,threshold, &
  !          & mat_c_columns, mat_c_rows, matAB)
  !   else
  !     call PruneList(blocked_memory_pool,alpha,threshold, &
  !          & mat_c_columns, mat_c_rows, matAB)
  !   end if
  !   !call StopTimer("Prune List")
  !
  !   if (present(beta_in) .and. abs(beta_in) .gt. 0) then
  !     call ScaleSparseMatrix(matC,beta)
  !     call IncrementSparseMatrix(matAB,matC)
  !   else
  !     call CopySparseMatrix(matAB,matC)
  !   end if
  !
  !   call DestructSparseMatrix(matAB)
  !   call DestructMatrixMemoryPool(blocked_memory_pool)
  !   !call StopTimer("Prune List")
  ! end subroutine Gemm2
  SUBROUTINE Gemm(matA, matB, matC, IsATransposed_in, IsBTransposed_in, &
       & alpha_in, beta_in, threshold_in, blocked_memory_pool_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in)  :: matA
    TYPE(SparseMatrix_t), INTENT(in)  :: matB
    TYPE(SparseMatrix_t), INTENT(inout) :: matC
    LOGICAL, OPTIONAL, INTENT(in) :: IsATransposed_in
    LOGICAL, OPTIONAL, INTENT(in) :: IsBTransposed_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    TYPE(MatrixMemoryPool_t), OPTIONAL, &
         & INTENT(inout), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(SparseMatrix_t) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(SparseMatrix_t) :: matAT, matBT
    TYPE(MatrixMemoryPool_t) :: blocked_memory_pool
    !! Counters and temporary data
    INTEGER :: mat_c_columns, mat_c_rows
    !! For Efficiency Purposes
    REAL(NTREAL) :: sparsity_estimate
    LOGICAL :: pool_flag

    !! Process Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(beta_in)) THEN
       beta = 0.0
    ELSE
       beta = beta_in
    END IF
    IF (.NOT. PRESENT(IsATransposed_in)) THEN
       IsATransposed = .FALSE.
    ELSE
       IsATransposed = IsATransposed_in
    END IF
    IF (.NOT. PRESENT(IsBTransposed_in)) THEN
       IsBTransposed = .FALSE.
    ELSE
       IsBTransposed = IsBTransposed_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Storage details for result matrix
    IF (IsATransposed) THEN
       mat_c_rows = matA%columns
    ELSE
       mat_c_rows = matA%rows
    END IF
    IF (IsBTransposed) THEN
       mat_c_columns = matB%rows
    ELSE
       mat_c_columns = matB%columns
    END IF

    !! Initialization of Memory
    !call StartTimer("Initialize Prune")
    !sparsity_estimate = max(dble(size(matA%values))/(matA%rows*matA%columns), &
    !  & dble(size(matB%values))/(matB%rows*matB%columns))
    sparsity_estimate = 4*MAX(DBLE(SIZE(matA%values))/(matA%rows*matA%columns),&
         & DBLE(SIZE(matB%values))/(matB%rows*matB%columns))
    IF (sparsity_estimate > 1.0) THEN
       sparsity_estimate = 1.0
    ELSE IF (sparsity_estimate < 1e-8) THEN
       sparsity_estimate = 1e-8
    END IF

    IF (.NOT. PRESENT(blocked_memory_pool_in)) THEN
       CALL ConstructMatrixMemoryPool(blocked_memory_pool,mat_c_columns, &
            & mat_c_rows, sparsity_estimate)
       pool_flag = .FALSE.
    ELSEIF (.NOT. CheckMemoryPoolValidity(blocked_memory_pool_in, &
         & mat_c_columns, mat_c_rows)) THEN
       CALL DestructMatrixMemoryPool(blocked_memory_pool_in)
       CALL ConstructMatrixMemoryPool(blocked_memory_pool_in,mat_c_columns, &
            & mat_c_rows, sparsity_estimate)
       pool_flag = .TRUE.
    ELSE
       CALL SetPoolSparsity(blocked_memory_pool_in, sparsity_estimate)
       pool_flag = .TRUE.
    END IF

    !! Block A and B
    !call StartTimer("Transposes")
    IF (.NOT. IsATransposed) THEN
       CALL TransposeSparseMatrix(matA,matAT)
    END IF
    IF (.NOT. IsBTransposed) THEN
       CALL TransposeSparseMatrix(matB,matBT)
    END IF
    !call StopTimer("Transposes")

    !call StopTimer("Initialize Prune")

    !call StartTimer("MM Loop")
    !call StartTimer("Multiply Block")
    IF (pool_flag) THEN
       IF (IsATransposed .AND. IsBTransposed) THEN
          CALL MultiplyBlock2(matA, matB, blocked_memory_pool_in)
       ELSEIF (IsATransposed) THEN
          CALL MultiplyBlock2(matA, matBT, blocked_memory_pool_in)
       ELSEIF (IsBTransposed) THEN
          CALL MultiplyBlock2(matAT, matB, blocked_memory_pool_in)
       ELSE
          CALL MultiplyBlock2(matAT, matBT, blocked_memory_pool_in)
       END IF
    ELSE
       IF (IsATransposed .AND. IsBTransposed) THEN
          CALL MultiplyBlock2(matA, matB, blocked_memory_pool)
       ELSEIF (IsATransposed) THEN
          CALL MultiplyBlock2(matA, matBT, blocked_memory_pool)
       ELSEIF (IsBTransposed) THEN
          CALL MultiplyBlock2(matAT, matB, blocked_memory_pool)
       ELSE
          CALL MultiplyBlock2(matAT, matBT, blocked_memory_pool)
       END IF
    END IF
    !call StopTimer("Multiply Block")

    !! Go from triplets to return matrix
    !call StartTimer("Prune List")
    !call StartTimer("Prune List")
    IF (pool_flag) THEN
       CALL PruneList2(blocked_memory_pool_in,alpha,threshold, &
            & mat_c_columns, mat_c_rows, matAB)
    ELSE
       CALL PruneList2(blocked_memory_pool,alpha,threshold, &
            & mat_c_columns, mat_c_rows, matAB)
    END IF
    !call StopTimer("Prune List")

    IF (PRESENT(beta_in)) THEN
       IF (ABS(beta_in) .GT. 0) THEN
          CALL ScaleSparseMatrix(matC,beta)
          CALL IncrementSparseMatrix(matAB,matC)
       ELSE
          CALL CopySparseMatrix(matAB,matC)
       END IF
    ELSE
       CALL CopySparseMatrix(matAB,matC)
    END IF

    CALL DestructSparseMatrix(matAB)
    CALL DestructMatrixMemoryPool(blocked_memory_pool)
    !call StopTimer("Prune List")
  END SUBROUTINE Gemm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  !! @param[in] this the matrix to compute the norm of.
  !! @param[out] norm_per_column the norm value for each column in this matrix.
  PURE SUBROUTINE SparseMatrixNorm(this, norm_per_column)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, &
         & INTENT(out) :: norm_per_column
    !! Local Data
    INTEGER :: outer_counter, inner_counter
    INTEGER :: elements_per_inner
    REAL(NTREAL) :: temp_value

    !! Allocate Space For Result
    ALLOCATE(norm_per_column(this%columns))
    norm_per_column = 0

    !! Iterate Over Local Data
    DO outer_counter = 1, this%columns
       elements_per_inner = this%outer_index(outer_counter+1) - &
            & this%outer_index(outer_counter)
       DO inner_counter = 1, elements_per_inner
          temp_value = this%values(this%outer_index(outer_counter)+ &
               & inner_counter)
          norm_per_column(outer_counter) = norm_per_column(outer_counter) + &
               & ABS(temp_value)
       END DO
    END DO
  END SUBROUTINE SparseMatrixNorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  !! @param[in] this the matrix to sum
  !! @result sum_value the sum of the matrix elements
  PURE FUNCTION SparseMatrixGrandSum(this) RESULT(sum_value)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    REAL(NTREAL) :: sum_value

    sum_value = SUM(this%values)

  END FUNCTION SparseMatrixGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE TransposeSparseMatrix(this, matT)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in)  :: this
    TYPE(SparseMatrix_t), INTENT(inout) :: matT
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    !! Temporary variables and counters
    INTEGER :: counter
    INTEGER :: temp_int
    INTEGER :: temp_columns, temp_rows
    INTEGER :: size_of_this

    CALL DestructSparseMatrix(matT)

    temp_rows = this%columns
    temp_columns = this%rows

    !! Construct a transposed triplet list
    size_of_this = this%outer_index(this%columns+1)
    CALL ConstructTripletList(triplet_list,size_of_this)
    CALL MatrixToTripletList(this,triplet_list)
    DO counter = 1,triplet_list%CurrentSize
       temp_int = triplet_list%data(counter)%index_column
       triplet_list%data(counter)%index_column = &
            & triplet_list%data(counter)%index_row
       triplet_list%data(counter)%index_row = temp_int
    END DO

    !! Build a new matrix from the triplet list
    CALL SortTripletList(triplet_list,temp_columns,sorted_triplet_list, &
         & bubble_in=.FALSE.)
    CALL ConstructFromTripletList(matT, sorted_triplet_list,&
         & temp_rows,temp_columns)

    !! Cleanup
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE TransposeSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !! the first matrix are followed by the columns of the matrices in the list.
  !! @param[in] mat_list list of matrices to compose.
  !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
  PURE SUBROUTINE ComposeSparseMatrixColumns(mat_list, out_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:), INTENT(in) :: mat_list
    TYPE(SparseMatrix_t), INTENT(inout) :: out_matrix
    !! Local Variables
    INTEGER :: total_columns, total_values
    INTEGER :: inner_start, inner_length
    INTEGER :: outer_start, outer_length
    INTEGER :: outer_offset
    INTEGER :: counter
    INTEGER :: size_of_mat

    CALL DestructSparseMatrix(out_matrix)

    !! Figure Out The Sizes
    total_columns = 0
    total_values  = 0
    DO counter = LBOUND(mat_list,dim=1), UBOUND(mat_list,dim=1)
       total_columns = total_columns + mat_list(counter)%columns
       size_of_mat = mat_list(counter)%outer_index(mat_list(counter)%columns+1)
       total_values  = total_values + size_of_mat
    END DO

    !! Allocate The Space
    CALL ConstructEmptySparseMatrix(out_matrix,total_columns,mat_list(1)%rows)
    ALLOCATE(out_matrix%inner_index(total_values))
    ALLOCATE(out_matrix%values(total_values))

    !! Fill In The Values
    inner_start = 1
    outer_start = 1
    outer_offset = 0
    DO counter = LBOUND(mat_list,dim=1),UBOUND(mat_list,dim=1)
       !! Inner indices and values
       size_of_mat = mat_list(counter)%outer_index(mat_list(counter)%columns+1)
       inner_length = size_of_mat
       out_matrix%inner_index(inner_start:inner_start+inner_length-1) = &
            & mat_list(counter)%inner_index
       out_matrix%values(inner_start:inner_start+inner_length-1) = &
            & mat_list(counter)%values
       inner_start = inner_start + inner_length
       !! Outer Indices
       !outer_length = size(mat_list(counter)%outer_index)
       outer_length = mat_list(counter)%columns+1
       out_matrix%outer_index(outer_start:outer_start+outer_length-1) = &
            & mat_list(counter)%outer_index + outer_offset
       outer_start = outer_start + outer_length - 1
       outer_offset = out_matrix%outer_index(outer_start)
    END DO

  END SUBROUTINE ComposeSparseMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a matrix, and split into into small blocks.
  !! @param[in] this matrix to perform this operation on.
  !! @param[in] num_blocks number of blocks to split into.
  !! @param[out] split_list 1D array of blocks.
  !! @param[out] block_offsets_out the offsets used for splitting.
  PURE SUBROUTINE SplitSparseMatrixColumns(this, num_blocks, split_list, &
       & block_offsets_out)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    INTEGER, INTENT(in) :: num_blocks
    TYPE(SparseMatrix_t), DIMENSION(num_blocks), INTENT(out) :: split_list
    INTEGER, DIMENSION(num_blocks+1), INTENT(out), OPTIONAL :: block_offsets_out
    !! Local Data
    INTEGER, DIMENSION(num_blocks) :: block_sizes
    INTEGER, DIMENSION(num_blocks+1) :: block_offsets
    INTEGER :: split_divisor
    !! Counters
    INTEGER :: split_counter
    !! Temporary variables
    INTEGER :: loffset, lcolumns, linner_offset, total_values

    !! Handle trivial case
    IF (num_blocks .EQ. 1) THEN
       CALL CopySparseMatrix(this,split_list(1))
       block_offsets(1) = 1
       block_offsets(2) = this%columns+1
    ELSE
       !! Compute the sizes of each block
       split_divisor = this%columns/num_blocks
       block_sizes = split_divisor
       !! Handle an uneven split
       block_sizes(num_blocks) = this%columns - split_divisor*(num_blocks-1)
       !! Offsets
       block_offsets(1) = 1
       DO split_counter = 2, num_blocks+1
          block_offsets(split_counter) = block_offsets(split_counter-1) + &
               & block_sizes(split_counter-1)
       END DO

       !! Split up the columns
       DO split_counter = 1, num_blocks
          !! Temporary variables
          loffset = block_offsets(split_counter)
          lcolumns = block_sizes(split_counter)
          linner_offset = this%outer_index(loffset)+1
          !! Construct
          CALL ConstructEmptySparseMatrix(split_list(split_counter), &
               & columns=lcolumns, rows=this%rows)
          !! Copy Outer Index
          split_list(split_counter)%outer_index =        &
               & this%outer_index(loffset:loffset+lcolumns)
          split_list(split_counter)%outer_index =        &
               & split_list(split_counter)%outer_index -    &
               & split_list(split_counter)%outer_index(1)
          total_values = split_list(split_counter)%outer_index(lcolumns+1)
          !! Copy Inner Indices and Values
          ALLOCATE(split_list(split_counter)%inner_index(total_values))
          split_list(split_counter)%inner_index = &
               & this%inner_index(linner_offset:linner_offset+total_values-1)
          ALLOCATE(split_list(split_counter)%values(total_values))
          split_list(split_counter)%values = &
               & this%values(linner_offset:linner_offset+total_values-1)
       END DO
    END IF
    IF (PRESENT(block_offsets_out)) THEN
       block_offsets_out = block_offsets
    END IF
  END SUBROUTINE SplitSparseMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  !! @param[in] this the matrix to construct the triplet list from.
  !! @param[out] triplet_list the triplet list we created.
  PURE SUBROUTINE MatrixToTripletList(this, triplet_list)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    TYPE(TripletList_t), INTENT(inout) :: triplet_list
    !! Helper variables
    INTEGER :: outer_counter, inner_counter
    INTEGER :: elements_per_inner
    INTEGER :: total_counter
    TYPE(Triplet_t) :: temporary
    INTEGER :: size_of_this

    size_of_this = this%outer_index(this%columns+1)
    CALL ConstructTripletList(triplet_list,size_of_this)

    total_counter = 1
    DO outer_counter = 1, this%columns
       elements_per_inner = this%outer_index(outer_counter+1) - &
            & this%outer_index(outer_counter)
       DO inner_counter = 1, elements_per_inner
          temporary%index_column = outer_counter
          temporary%index_row = this%inner_index(total_counter)
          temporary%point_value = this%values(total_counter)
          triplet_list%data(total_counter) = temporary
          total_counter = total_counter + 1
       END DO
    END DO
  END SUBROUTINE MatrixToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! We first create a triplet list, and then call the print triplet list
  !! function.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintSparseMatrix(this, file_name_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: file_name_in
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    INTEGER :: file_handler
    INTEGER :: counter
    INTEGER :: size_of_this

    !! Process Optional Parameters
    IF (PRESENT(file_name_in)) THEN
       file_handler = 16
       OPEN(unit = file_handler, file = file_name_in)
    ELSE
       file_handler = 6
    END IF

    !! Print
    CALL MatrixToTripletList(this,triplet_list)

    size_of_this = this%outer_index(this%columns+1)

    WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate real general"
    WRITE(file_handler,'(A)') "%"
    WRITE(file_handler,*) this%rows, this%columns, size_of_this
    DO counter = 1,size_of_this
       WRITE(file_handler,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value
    END DO

    IF (PRESENT(file_name_in)) THEN
       CLOSE(file_handler)
    END IF
    CALL DestructTripletList(triplet_list)
  END SUBROUTINE PrintSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if a matrix is equal to the identity matrix.
  !! This routine is really just for testing. You can multiply a matrix and
  !! its inverse, and then call this routine to make sure multiplication is
  !! correct.
  !! @param[in] this the matrix to check.
  !! @param[in] threshold_in for flushing values to zero. Default value is 10^-8.
  !! @result true if the matrix is equal to the identity matrix.
  PURE FUNCTION CheckIfIdentity(this, threshold_in) RESULT(isidentity)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: this
    LOGICAL :: isidentity
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    INTEGER :: counter
    TYPE(Triplet_t) :: temp
    REAL(NTREAL) :: threshold

    !! Process Optional Parameters
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 1e-8
    ELSE
       threshold = threshold_in
    END IF

    CALL MatrixToTripletList(this,triplet_list)
    isidentity = .TRUE.
    DO counter=1, triplet_list%CurrentSize
       temp = triplet_list%data(counter)
       IF (temp%index_row .EQ. temp%index_column) THEN
          IF (ABS(temp%point_value - 1.0) > threshold) THEN
             isidentity = .FALSE.
             EXIT
          ENDIF
       ELSE
          IF (ABS(temp%point_value) > threshold) THEN
             isidentity = .FALSE.
             EXIT
          ENDIF
       ENDIF
    END DO
    CALL DestructTripletList(triplet_list)
  END FUNCTION CheckIfIdentity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! pure subroutine MultiplyBlock(matAT,matBT,value_array_ptr)
  !   !! Parameters
  !   type(SparseMatrix_t), intent(in)  :: matAT
  !   type(SparseMatrix_t), intent(in)  :: matBT
  !   real(NTREAL), dimension(:,:), intent(inout) :: value_array_ptr
  !   !! Temp Variables
  !   real(NTREAL) :: temp_value_a, temp_value_b
  !   integer :: temp_index_a, temp_index_b
  !   integer :: elements_per_inner_a
  !   integer :: elements_per_inner_b
  !   !! Counters
  !   integer :: outer_counter, inner_counter_a, inner_counter_b
  !
  !   !! Multiply
  !   do outer_counter = 1, matAT%columns
  !     elements_per_inner_a = matAT%outer_index(outer_counter+1) - &
  !       & matAT%outer_index(outer_counter)
  !     do inner_counter_a = 1, elements_per_inner_a
  !       temp_value_a = matAT%values(matAT%outer_index(outer_counter)+ &
  !         & inner_counter_a)
  !       temp_index_a = matAT%inner_index(matAT%outer_index(outer_counter)+ &
  !           & inner_counter_a)
  !       elements_per_inner_b = matBT%outer_index(temp_index_a+1) - &
  !         & matBT%outer_index(temp_index_a)
  !       do inner_counter_b = 1, elements_per_inner_b
  !         temp_index_b = matBT%inner_index(matBT%outer_index(temp_index_a)+ &
  !           & inner_counter_b)
  !         temp_value_b = matBT%values(matBT%outer_index(temp_index_a)+ &
  !           & inner_counter_b)
  !         value_array_ptr(temp_index_b,outer_counter) = &
  !           & value_array_ptr(temp_index_b,outer_counter) + &
  !           & temp_value_a*temp_value_b
  !       end do
  !     end do
  !   end do
  ! end subroutine MultiplyBlock
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE MultiplyBlock2(matAT,matBT,memorypool)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in)  :: matAT
    TYPE(SparseMatrix_t), INTENT(in)  :: matBT
    TYPE(MatrixMemoryPool_t), INTENT(inout) :: memorypool
    !! Temp Variables
    REAL(NTREAL) :: temp_value_a, temp_value_b, temp_value_c
    INTEGER :: temp_inserted_values
    INTEGER :: temp_index_a, temp_index_b
    INTEGER :: elements_per_inner_a
    INTEGER :: elements_per_inner_b
    !! Counters
    INTEGER :: outer_counter, inner_counter_a, inner_counter_b

    !! Multiply
    DO outer_counter = 1, matAT%columns
       elements_per_inner_a = matAT%outer_index(outer_counter+1) - &
            & matAT%outer_index(outer_counter)
       DO inner_counter_a = 1, elements_per_inner_a
          temp_value_a = matAT%values(matAT%outer_index(outer_counter)+ &
               & inner_counter_a)
          temp_index_a = matAT%inner_index(matAT%outer_index(outer_counter)+ &
               & inner_counter_a)
          elements_per_inner_b = matBT%outer_index(temp_index_a+1) - &
               & matBT%outer_index(temp_index_a)
          DO inner_counter_b = 1, elements_per_inner_b
             temp_index_b = matBT%inner_index(matBT%outer_index(temp_index_a)+ &
                  & inner_counter_b)
             temp_value_b = matBT%values(matBT%outer_index(temp_index_a)+ &
                  & inner_counter_b)
             temp_value_c = memorypool%value_array(temp_index_b,outer_counter)
             IF (temp_value_c .EQ. 0) THEN
                ! temp_inserted_values = memorypool%inserted_per_bucket_ptr(&
                !   & temp_index_b/memorypool%hash_size+1,outer_counter) + 1
                temp_inserted_values = memorypool%inserted_per_bucket(&
                     & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) + 1
                ! memorypool%inserted_per_bucket_ptr(&
                !   & temp_index_b/memorypool%hash_size+1,outer_counter) = &
                !   & temp_inserted_values
                memorypool%inserted_per_bucket(&
                     & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) = &
                     & temp_inserted_values
                ! memorypool%hash_index_ptr(temp_inserted_values, &
                !   & temp_index_b/memorypool%hash_size+1,outer_counter) = &
                !   & temp_index_b
                memorypool%hash_index(temp_inserted_values+ &
                     & ((temp_index_b-1)/memorypool%hash_size)*memorypool%hash_size, &
                     & outer_counter) = temp_index_b
             END IF
             memorypool%value_array(temp_index_b,outer_counter) = &
                  & temp_value_c + temp_value_a*temp_value_b
          END DO
       END DO
    END DO
  END SUBROUTINE MultiplyBlock2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! pure subroutine PruneList(blocked_memory_pool,alpha,threshold, &
  !      & mat_c_columns, mat_c_rows, matAB)
  !   !! Parameters
  !   type(MatrixMemoryPool_t), intent(inout) :: blocked_memory_pool
  !   real(NTREAL), intent(in) :: alpha
  !   real(NTREAL), intent(in) :: threshold
  !   integer, intent(in) :: mat_c_columns
  !   integer, intent(in) :: mat_c_rows
  !   type(SparseMatrix_t), intent(inout) :: matAB
  !   !! Local data
  !   integer :: pruned_counter
  !   integer :: row_counter_c, column_counter_c
  !   real(NTREAL) :: working_value
  !   type(TripletList_t) :: unsorted_pruned_list
  !   type(TripletList_t) :: sorted_pruned_list
  !
  !   pruned_counter = 1
  !   do row_counter_c = 1, mat_c_rows
  !     do column_counter_c = 1, mat_c_columns
  !       working_value = blocked_memory_pool%value_array(column_counter_c,&
  !         & row_counter_c)
  !       blocked_memory_pool%value_array(column_counter_c, &
  !         & row_counter_c) = 0
  !       if (abs(working_value) .gt. threshold) then
  !         blocked_memory_pool%pruned_list(pruned_counter)%point_value = &
  !           & alpha*working_value
  !         blocked_memory_pool%pruned_list(pruned_counter)%index_column = &
  !           & column_counter_c
  !         blocked_memory_pool%pruned_list(pruned_counter)%index_row = &
  !           & row_counter_c
  !         pruned_counter = pruned_counter + 1
  !       end if
  !     end do
  !   end do
  !   call ConstructTripletList(unsorted_pruned_list,pruned_counter-1)
  !   unsorted_pruned_list%data = &
  !     & blocked_memory_pool%pruned_list(1:pruned_counter-1)
  !   call SortTripletList(unsorted_pruned_list,mat_c_columns,sorted_pruned_list)
  !   !call SortTripletList(blocked_memory_pool%pruned_list(1:pruned_counter-1),&
  !   !     & mat_c_columns, sorted_pruned_list)
  !   call ConstructFromTripletList(matAB,sorted_pruned_list,mat_c_rows, &
  !        & mat_c_columns)
  !   call DestructTripletList(sorted_pruned_list)
  ! end subroutine PruneList
  PURE SUBROUTINE PruneList2(memorypool,alpha,threshold, &
       & mat_c_columns, mat_c_rows, matAB)
    !! Parameters
    TYPE(MatrixMemoryPool_t), INTENT(inout) :: memorypool
    REAL(NTREAL), INTENT(in) :: alpha
    REAL(NTREAL), INTENT(in) :: threshold
    INTEGER, INTENT(in) :: mat_c_columns
    INTEGER, INTENT(in) :: mat_c_rows
    TYPE(SparseMatrix_t), INTENT(inout) :: matAB
    !! Local data
    INTEGER :: row_counter_c, column_counter_c, hash_counter
    REAL(NTREAL) :: working_value
    INTEGER :: working_column
    INTEGER :: temp_values_per_hash
    TYPE(TripletList_t) :: unsorted_pruned_list
    TYPE(TripletList_t) :: sorted_pruned_list
    INTEGER :: pruned_counter

    pruned_counter = 1
    DO row_counter_c = 1, mat_c_rows
       DO column_counter_c = 1, (mat_c_columns-1)/memorypool%hash_size+1
          !! Sort the elements in a hash
          ! temp_values_per_hash = memorypool%inserted_per_bucket_ptr(&
          !   & column_counter_c, row_counter_c)
          temp_values_per_hash = memorypool%inserted_per_bucket(&
               & column_counter_c,row_counter_c)
          ! memorypool%inserted_per_bucket_ptr(column_counter_c, row_counter_c) = 0
          memorypool%inserted_per_bucket(column_counter_c,row_counter_c) = 0
          !! Copy them
          DO hash_counter=1,temp_values_per_hash
             ! working_column = memorypool%hash_index_ptr(hash_counter, &
             !   & column_counter_c,row_counter_c)
             working_column = memorypool%hash_index(hash_counter+ &
                  & (column_counter_c-1)*memorypool%hash_size, row_counter_c)
             working_value = memorypool%value_array(working_column,row_counter_c)
             memorypool%value_array(working_column,row_counter_c) = 0
             IF (ABS(alpha*working_value) .GT. threshold) THEN
                memorypool%pruned_list(pruned_counter)%point_value = &
                     & alpha*working_value
                memorypool%pruned_list(pruned_counter)%index_column = &
                     & working_column
                memorypool%pruned_list(pruned_counter)%index_row = &
                     & row_counter_c
                pruned_counter = pruned_counter + 1
             END IF
          END DO
       END DO
    END DO
    CALL ConstructTripletList(unsorted_pruned_list,pruned_counter-1)
    unsorted_pruned_list%data = memorypool%pruned_list(1:pruned_counter-1)
    CALL SortTripletList(unsorted_pruned_list,mat_c_columns,sorted_pruned_list,&
         & bubble_in=.TRUE.)
    CALL ConstructFromTripletList(matAB,sorted_pruned_list,mat_c_rows, &
         & mat_c_columns)
    CALL DestructTripletList(sorted_pruned_list)
  END SUBROUTINE PruneList2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SparseMatrixModule
