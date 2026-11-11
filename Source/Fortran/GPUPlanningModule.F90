! Pure, hardware-independent estimates; runtime memory admission is separate.
MODULE GPUPlanningModule
  USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EstimateGPUPlan, GPUMemoryBudget
CONTAINS
  PURE SUBROUTINE EstimateGPUPlan(n,pr,pc,nnz_wh,nnz_x,iterations,peak_bytes,seconds,index_safe)
    INTEGER, INTENT(IN) :: n,pr,pc,iterations
    INTEGER(INT64), INTENT(IN) :: nnz_wh,nnz_x
    INTEGER(INT64), INTENT(OUT) :: peak_bytes
    REAL(REAL64), INTENT(OUT) :: seconds
    LOGICAL, INTENT(OUT) :: index_safe
    REAL(REAL64) :: rows,columns,cells,peers,degree,products,peak,p,runtime
    peak_bytes=HUGE(peak_bytes); seconds=HUGE(seconds); index_safe=.FALSE.
    IF(n<1 .OR. pr<1 .OR. pc<1 .OR. iterations<1 .OR. nnz_wh<0 .OR. nnz_x<0) RETURN
    IF(MOD(n,pr)/=0 .OR. MOD(n,pc)/=0) RETURN
    rows=REAL(n/pr,REAL64); columns=REAL(n/pc,REAL64); cells=rows*columns
    p=REAL(pr,REAL64)*REAL(pc,REAL64)
    IF(p>REAL(HUGE(0),REAL64)) RETURN
    ! All current CSR pointers/counts are 32-bit, including intermediate panels.
    ! Reject layouts whose worst-case dense result or gathered panels overflow.
    IF(MAX(cells,cells*pr,cells*pc)>REAL(HUGE(0),REAL64)) RETURN
    ! GPU X, X^2 and merge scratch: at most three full local CSR blocks.
    ! Panel payloads: local + row receive + column receive + assembled A.
    peers=1.0_REAL64+2.0_REAL64*pc+pr
    peak=(36.0_REAL64+12.0_REAL64*peers+8.0_REAL64)*cells
    peak=peak+4.0_REAL64*rows*CEILING(columns/32.0_REAL64,KIND=INT64)
    ! WH cannot fill in during the region; use its global nnz as a safe bound
    ! on every rank's local storage, rather than assume balanced input.
    peak=peak+12.0_REAL64*REAL(nnz_wh,REAL64)
    ! Row pointers, scans/reductions, allocator/protocol workspace and padding.
    peak=peak+512.0_REAL64*(REAL(n,REAL64)+rows+columns+1.0_REAL64)+1073741824.0_REAL64
    IF(peak>=REAL(HUGE(peak_bytes),REAL64)) RETURN
    peak_bytes=CEILING(peak,KIND=INT64); index_safe=.TRUE.
    ! Ranking model for a native API call, including communicator setup.
    ! Coefficients are initial fits to the measured 10k/64k regimes, not a
    ! promised runtime. Degree-based work is bounded by a dense row.
    degree=MIN(REAL(n,REAL64),REAL(nnz_x,REAL64)/REAL(n,REAL64))
    products=REAL(n,REAL64)*degree*degree
    runtime=0.0007_REAL64+2e-12_REAL64*cells+5e-12_REAL64*products/p
    runtime=runtime+1e-10_REAL64*REAL(nnz_x,REAL64)*(REAL(pr,REAL64)+REAL(pc,REAL64))/p
    seconds=REAL(iterations,REAL64)*runtime+0.02_REAL64
    IF(p>1.0_REAL64) THEN
      seconds=seconds+2.0_REAL64+REAL(iterations,REAL64)*1e-4_REAL64*(REAL(pr,REAL64)+REAL(pc,REAL64))
    ELSE
      seconds=seconds+24.0_REAL64*(REAL(nnz_wh,REAL64)+REAL(nnz_x,REAL64))/1e10_REAL64
    END IF
  END SUBROUTINE

  PURE FUNCTION GPUMemoryBudget(free_bytes,total_bytes,limit_bytes) RESULT(budget)
    INTEGER(INT64), INTENT(IN) :: free_bytes,total_bytes,limit_bytes
    INTEGER(INT64) :: budget
    ! 20% of currently free memory is left for runtime and external users.
    ! A user cap can reduce admission; it never increases live free capacity.
    budget=INT(0.8_REAL64*REAL(MAX(0_INT64,MIN(free_bytes,total_bytes)),REAL64),INT64)
    IF(limit_bytes>=0) budget=MIN(budget,limit_bytes)
  END FUNCTION
END MODULE
