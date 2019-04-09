!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains some enumerators which name the tasks for Gemm.
MODULE GemmTasksModule
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(c)
    !> Something is in progress
    ENUMERATOR :: TaskRunningA
    !> First we gather the blocks of A and send the size.
    ENUMERATOR :: LocalGatherA
    !> After the local gather, we then send the size to the other tasks.
    ENUMERATOR :: SendSizeA
    !> Next we compose those blocks of A into one big send buffer and send.
    ENUMERATOR :: ComposeA
    !> Wait for the outer index values to be gathered.
    ENUMERATOR :: WaitOuterA
    !> Wait for the inner index values to be gathered.
    ENUMERATOR :: WaitInnerA
    !> Wait for the data values to be gathered,
    ENUMERATOR :: WaitDataA
    !> Need to adjusts indices, transpose the values of A.
    ENUMERATOR :: AdjustIndicesA
    !> Just waiting on that last task.
    ENUMERATOR :: CleanupA
    !> No more work to do.
    ENUMERATOR :: FinishedA
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(c)
    !> Something is in progress
    ENUMERATOR :: TaskRunningB
    !> First we gather the blocks of B and send the size.
    ENUMERATOR :: LocalGatherB
    !> Next we compose those blocks of B into one big send buffer and send.
    ENUMERATOR :: LocalComposeB
    !> After the local gather, we then send the size to the other tasks.
    ENUMERATOR :: SendSizeB
    !> Wait for the outer index values to be gathered.
    ENUMERATOR :: WaitOuterB
    !> Wait for the inner index values to be gathered.
    ENUMERATOR :: WaitInnerB
    !> Wait for the data values to be gathered, and then adjusts the indices.
    ENUMERATOR :: WaitDataB
    !> Need to adjusts indices of B.
    ENUMERATOR :: AdjustIndicesB
    !> Just waiting on that last task.
    ENUMERATOR :: CleanupB
    !> No more work to do.
    ENUMERATOR :: FinishedB
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(c)
    !> Something is in progress.
    ENUMERATOR :: TaskRunningAB
    !> A and B matrix both missing, so it cannot do gemm.
    ENUMERATOR :: AwaitingAB
    !> Actually call gemm and compute a block, and send its size.
    ENUMERATOR :: GemmAB
    !> After the local Gemm, we then send the size to the other tasks.
    ENUMERATOR :: SendSizeAB
    !> Start sending the data for summing.
    ENUMERATOR :: GatherAndSumAB
    !> Wait for the outer index values to be gathered.
    ENUMERATOR :: WaitOuterAB
    !> Wait for the inner index values to be gathered.
    ENUMERATOR :: WaitInnerAB
    !> Wait for the data values to be gathered. Once receive, we increment.
    ENUMERATOR :: WaitDataAB
    !> Sum up the gathered matrices.
    ENUMERATOR :: LocalSumAB
    !> Just waiting on that last task.
    ENUMERATOR :: CleanupAB
    !> No more work to do.
    enumerator :: FinishedAB
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE GemmTasksModule
