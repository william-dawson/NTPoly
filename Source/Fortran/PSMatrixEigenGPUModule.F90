! Optional real FP64 cuSOLVERMp boundary for the ordinary dense eigen API.
MODULE PSMatrixEigenGPUModule
  USE ISO_C_BINDING
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE DataTypesModule, ONLY: NTREAL
  USE NTMPIModule
  USE PSMatrixModule, ONLY: Matrix_ps, GetMatrixTripletList, ConstructEmptyMatrix, FillMatrixFromTripletList
  USE TripletListModule, ONLY: TripletList_r, ConstructTripletList, DestructTripletList, AppendToTripletList
  USE TripletModule, ONLY: Triplet_r
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TryGPUEigenDecomposition
  INTERFACE
    INTEGER(C_INT) FUNCTION solve(comm,n,pr,pc,count,rows,columns,values,nvals,vectors,threshold,root,output) &
        BIND(C,NAME='ntpoly_cusolvermp_eigh_')
      IMPORT C_INT,C_DOUBLE,C_PTR
      INTEGER(C_INT) :: comm,n,pr,pc,count,rows(*),columns(*),nvals,vectors,root
      REAL(C_DOUBLE) :: values(*),threshold
      TYPE(C_PTR) :: output
    END FUNCTION
    SUBROUTINE result_info(output,count,values,eigen_count) BIND(C,NAME='ntpoly_cusolvermp_result_')
      IMPORT C_INT,C_PTR
      TYPE(C_PTR) :: output,values
      INTEGER(C_INT) :: count,eigen_count
    END SUBROUTINE
    SUBROUTINE coordinates(output,rows,columns,values) BIND(C,NAME='ntpoly_cusolvermp_coordinates_')
      IMPORT C_INT,C_DOUBLE,C_PTR
      TYPE(C_PTR) :: output
      INTEGER(C_INT) :: rows(*),columns(*)
      REAL(C_DOUBLE) :: values(*)
    END SUBROUTINE
    SUBROUTINE free_result(output) BIND(C,NAME='ntpoly_cusolvermp_free_')
      IMPORT C_PTR
      TYPE(C_PTR) :: output
    END SUBROUTINE
  END INTERFACE
CONTAINS
  SUBROUTINE TryGPUEigenDecomposition(matrix,eigenvalues,nvals,threshold,performed,eigenvectors)
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    TYPE(Matrix_ps), INTENT(INOUT) :: eigenvalues
    INTEGER, INTENT(IN) :: nvals
    REAL(NTREAL), INTENT(IN) :: threshold
    LOGICAL, INTENT(OUT) :: performed
    TYPE(Matrix_ps), OPTIONAL, INTENT(INOUT) :: eigenvectors
    TYPE(TripletList_r) :: entries,output_entries
    TYPE(Triplet_r) :: entry
    INTEGER(C_INT), ALLOCATABLE :: rows(:),columns(:)
    REAL(C_DOUBLE), ALLOCATABLE :: values(:)
    REAL(C_DOUBLE), POINTER :: ev(:)
    TYPE(C_PTR) :: output,ev_pointer
    CHARACTER(32) :: option,policy
    INTEGER :: flags(4),minimum(4),maximum(4),ierr,i
    INTEGER(C_INT) :: n,pr,pc,count,wanted,vectors,root,status,eigen_count
    performed=.FALSE.
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_EIGEN',option)
    flags=0
    IF(TRIM(option)=='1') flags(1)=1
    CALL GET_ENVIRONMENT_VARIABLE('NTPOLY_GPU_EIGEN_RANK_POLICY',policy)
    SELECT CASE(TRIM(policy))
    CASE('','auto'); flags(2)=0
    CASE('all'); flags(2)=1
    CASE('root'); flags(2)=2
    CASE DEFAULT; flags(2)=-1
    END SELECT
    flags(3)=nvals
    IF(PRESENT(eigenvectors)) flags(4)=1
    CALL MPI_Allreduce(flags,minimum,4,MPI_INTEGER,MPI_MIN,matrix%process_grid%global_comm,ierr)
    CALL MPI_Allreduce(flags,maximum,4,MPI_INTEGER,MPI_MAX,matrix%process_grid%global_comm,ierr)
    IF(ANY(minimum/=maximum)) THEN
      PRINT *, 'NTPoly GPU eigen: inconsistent request across ranks'
      CALL MPI_Abort(matrix%process_grid%global_comm,1,ierr)
    END IF
    IF(flags(1)==0) RETURN
    IF(flags(2)<0) THEN
      PRINT *, 'NTPoly GPU eigen: rank policy must be auto, all or root'
      CALL MPI_Abort(matrix%process_grid%global_comm,1,ierr)
    END IF
    IF(matrix%is_complex .OR. matrix%process_grid%num_process_slices/=1) RETURN
    IF(.NOT.IEEE_IS_FINITE(threshold) .OR. threshold<0 .OR. nvals<0) RETURN
    n=matrix%actual_matrix_dimension
    IF(n<1) RETURN
    wanted=MIN(nvals,n);vectors=flags(4)
    root=0
    IF(flags(2)==2 .OR. (flags(2)==0 .AND. n<=16384)) root=1
    IF(wanted==0) THEN
      CALL ConstructEmptyMatrix(eigenvalues,matrix)
      IF(PRESENT(eigenvectors)) CALL ConstructEmptyMatrix(eigenvectors,matrix)
      performed=.TRUE.;RETURN
    END IF
    CALL GetMatrixTripletList(matrix,entries)
    count=entries%CurrentSize
    ALLOCATE(rows(MAX(1,count)),columns(MAX(1,count)),values(MAX(1,count)))
    DO i=1,count
      rows(i)=entries%DATA(i)%index_row;columns(i)=entries%DATA(i)%index_column;values(i)=entries%DATA(i)%point_value
    END DO
    CALL DestructTripletList(entries)
    pr=matrix%process_grid%num_process_rows;pc=matrix%process_grid%num_process_columns
    status=solve(matrix%process_grid%global_comm,n,pr,pc,count,rows,columns,values,wanted,vectors,threshold,root,output)
    IF(status/=0) THEN
      IF(matrix%process_grid%global_rank==0) PRINT *, 'NTPoly cuSOLVERMp eigen: CPU fallback (capacity/layout unavailable)'
      RETURN
    END IF
    CALL result_info(output,count,ev_pointer,eigen_count)
    CALL ConstructTripletList(output_entries)
    IF(eigen_count>0) THEN
      CALL C_F_POINTER(ev_pointer,ev,[eigen_count])
      DO i=1,eigen_count
        IF(ABS(ev(i))<=threshold) CYCLE
        entry%index_row=i;entry%index_column=i;entry%point_value=ev(i)
        CALL AppendToTripletList(output_entries,entry)
      END DO
    END IF
    CALL ConstructEmptyMatrix(eigenvalues,matrix)
    CALL FillMatrixFromTripletList(eigenvalues,output_entries,preduplicated_in=.TRUE.)
    CALL DestructTripletList(output_entries)
    IF(PRESENT(eigenvectors)) THEN
      DEALLOCATE(rows,columns,values)
      ALLOCATE(rows(MAX(1,count)),columns(MAX(1,count)),values(MAX(1,count)))
      CALL coordinates(output,rows,columns,values)
      CALL ConstructTripletList(output_entries,count)
      DO i=1,count
        output_entries%DATA(i)%index_row=rows(i);output_entries%DATA(i)%index_column=columns(i)
        output_entries%DATA(i)%point_value=values(i)
      END DO
      CALL ConstructEmptyMatrix(eigenvectors,matrix)
      CALL FillMatrixFromTripletList(eigenvectors,output_entries,preduplicated_in=.TRUE.)
      CALL DestructTripletList(output_entries)
    END IF
    CALL free_result(output)
    performed=.TRUE.
  END SUBROUTINE
END MODULE
