PROGRAM TestGPUPlanning
  USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE GPUPlanningModule
  IMPLICIT NONE
  INTEGER(INT64) :: root,full,big,invalid
  REAL(REAL64) :: sparse_time,dense_time,full_time,t
  LOGICAL :: safe
  CALL EstimateGPUPlan(10000,1,1,400000_INT64,400000_INT64,24,root,sparse_time,safe)
  IF(.NOT.safe .OR. root<8_INT64*10000*10000) STOP 1
  CALL EstimateGPUPlan(10000,2,2,400000_INT64,400000_INT64,24,full,full_time,safe)
  IF(.NOT.safe .OR. full>=root .OR. sparse_time>=full_time) STOP 2
  ! Same dimension, high fill: work rather than a size cutoff changes the ranking.
  CALL EstimateGPUPlan(10000,1,1,100000000_INT64,100000000_INT64,24,big,dense_time,safe)
  CALL EstimateGPUPlan(10000,2,2,100000000_INT64,100000000_INT64,24,big,full_time,safe)
  IF(.NOT.safe .OR. dense_time<=full_time) STOP 3
  ! A cap between the layouts admits distributed execution, not coalescing.
  IF(GPUMemoryBudget(1000000000000_INT64,1000000000000_INT64,6000_INT64*1048576)<full) STOP 4
  IF(GPUMemoryBudget(1000000000000_INT64,1000000000000_INT64,6000_INT64*1048576)>=root) STOP 5
  IF(GPUMemoryBudget(100_INT64,100_INT64,-1_INT64)/=80) STOP 6
  IF(GPUMemoryBudget(100_INT64,100_INT64,200_INT64)/=80) STOP 7
  IF(GPUMemoryBudget(0_INT64,100_INT64,-1_INT64)/=0) STOP 8
  CALL EstimateGPUPlan(64000,2,2,2600000_INT64,2600000_INT64,24,big,t,safe)
  IF(.NOT.safe .OR. big<=full .OR. .NOT.IEEE_IS_FINITE(t)) STOP 9
  ! Boundaries reject malformed layouts and 32-bit panel/result overflow.
  CALL EstimateGPUPlan(64000,1,1,2600000_INT64,2600000_INT64,24,invalid,t,safe)
  IF(safe) STOP 10
  CALL EstimateGPUPlan(HUGE(0),1,1,HUGE(0_INT64),HUGE(0_INT64),24,invalid,t,safe)
  IF(safe .OR. invalid/=HUGE(invalid)) STOP 11
  CALL EstimateGPUPlan(10001,2,2,1_INT64,1_INT64,24,invalid,t,safe)
  IF(safe) STOP 12
  CALL EstimateGPUPlan(10000,0,2,1_INT64,1_INT64,24,invalid,t,safe)
  IF(safe) STOP 13
  PRINT *, 'GPU planning: fill-sensitive ranking, peak bounds, caps and integer overflow: PASS'
END PROGRAM
