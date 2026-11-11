PROGRAM TestSolverCSR
  USE SolverCSR
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: rp(:), ci(:), tp(:), ti(:), ends(:), offsets(:)
  INTEGER, ALLOCATABLE :: packed_rp(:), packed_ci(:)
  INTEGER, ALLOCATABLE :: partial_ci(:)
  REAL(C_DOUBLE), ALLOCATABLE :: partial_v(:)
  REAL(C_DOUBLE), ALLOCATABLE :: val(:), tv(:), packed_v(:)
  REAL(C_DOUBLE) :: dense(3,5), rebuilt(3,5)
  INTEGER :: r, j, ch, b, lo, at
  LOGICAL :: ok
  ! Exceed the production seed capacity with no initialized entries, then
  ! grow again while preserving a small live prefix and ignoring the tail.
  ALLOCATE(partial_ci(65536), partial_v(65536))
  CALL ReservePartials(0, 800000, partial_ci, partial_v)
  IF (SIZE(partial_ci) < 800000) ERROR STOP 'initial partial capacity'
  partial_ci(1:3) = [2,4,6]
  partial_v(1:3) = [1D0,2D0,3D0]
  CALL ReservePartials(3, 1600001, partial_ci, partial_v)
  IF (ANY(partial_ci(1:3)/=[2,4,6]) .OR. &
       ANY(partial_v(1:3)/=[1D0,2D0,3D0])) ERROR STOP 'partial growth prefix'
  rp = [0,3,3,6]
  ci = [1,2,5,2,3,5]
  val = [1D0,2D0,3D0,4D0,5D0,6D0]
  dense = 0
  DO r=1,3
    DO j=rp(r)+1,rp(r+1)
      dense(r,ci(j))=val(j)
    END DO
  END DO
  CALL TransposeCSR(3,5,rp,ci,val,tp,ti,tv)
  CALL TransposeCSR(5,3,tp,ti,tv,packed_rp,packed_ci,packed_v)
  IF (ANY(packed_rp/=rp) .OR. ANY(packed_ci/=ci) .OR. ANY(packed_v/=val)) ERROR STOP 'transpose round trip'
  DO b=1,2
    IF (b==1) THEN
      ends=[0,5]
    ELSE
      ends=[0,1,3,4,5] ! ragged slabs, empty slab, empty matrix row
    END IF
    CALL PackChunks(3,rp,ci,val,ends,packed_rp,packed_ci,packed_v,offsets)
    rebuilt=0
    DO ch=1,SIZE(ends)-1
      lo=(ch-1)*4
      DO r=1,3
        DO j=packed_rp(lo+r)+1,packed_rp(lo+r+1)
          at=offsets(ch)+j
          rebuilt(r,packed_ci(at)+ends(ch))=packed_v(at)
        END DO
      END DO
    END DO
    IF (ANY(rebuilt/=dense)) ERROR STOP 'chunk packing'
  END DO
  CALL PlanChunks(5,ci,[0,2,3,4,4,6],4_C_INT64_T,ends,ok)
  IF (.NOT.ok .OR. ANY(ends/=[0,2,4,5])) ERROR STOP 'workspace slabs'
  CALL PlanChunks(5,ci,[0,2,3,4,4,6],3_C_INT64_T,ends,ok)
  IF (ok) ERROR STOP 'single-k budget rejection'
  rp=0*rp
  ci=[INTEGER::]
  val=[REAL(C_DOUBLE)::]
  CALL TransposeCSR(3,5,rp,ci,val,tp,ti,tv)
  IF (ANY(tp/=0) .OR. SIZE(ti)/=0) ERROR STOP 'empty transpose'
  CALL PlanChunks(5,ci,[0,0,0,0,0,0],1_C_INT64_T,ends,ok)
  IF (.NOT.ok .OR. ANY(ends/=[0,5])) ERROR STOP 'empty chunks'
  CALL PackChunks(3,rp,ci,val,ends,packed_rp,packed_ci,packed_v,offsets)
  IF (SIZE(packed_ci)/=0 .OR. ANY(packed_rp/=0)) ERROR STOP 'empty packing'
  PRINT *, 'CSR transpose, ragged/empty packing, and workspace planning: PASS'
END PROGRAM
