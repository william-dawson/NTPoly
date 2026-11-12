! Host-side CSR transformations used by the purification driver.
! Pointers are zero-based; column indices are one-based (NTPoly convention).
MODULE SolverCSR
  USE ISO_C_BINDING, ONLY: C_DOUBLE, C_INT64_T
  IMPLICIT NONE
CONTAINS
  ! Preserve only the initialized prefix. In particular, the first product
  ! must not copy the unused/uninitialized contents of the seed allocation.
  SUBROUTINE ReservePartials(used, needed, ci, val)
    INTEGER, INTENT(IN) :: used, needed
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: ci(:)
    REAL(C_DOUBLE), ALLOCATABLE, INTENT(INOUT) :: val(:)
    INTEGER, ALLOCATABLE :: grown_ci(:)
    REAL(C_DOUBLE), ALLOCATABLE :: grown_val(:)
    INTEGER :: capacity
    IF (used < 0 .OR. needed < used) ERROR STOP 'ReservePartials: invalid size'
    capacity = 0
    IF (ALLOCATED(ci)) capacity = SIZE(ci)
    IF (needed <= capacity) RETURN
    IF (capacity > HUGE(capacity)/2) THEN
      capacity = needed
    ELSE
      capacity = MAX(needed, 2*capacity)
    END IF
    ALLOCATE(grown_ci(capacity), grown_val(capacity))
    IF (used > 0) THEN
      grown_ci(1:used) = ci(1:used)
      grown_val(1:used) = val(1:used)
    END IF
    CALL MOVE_ALLOC(grown_ci, ci)
    CALL MOVE_ALLOC(grown_val, val)
  END SUBROUTINE

  SUBROUTINE TransposeCSR(m, n, rp, ci, val, tp, ti, tv, phase_times)
    INTEGER, INTENT(IN) :: m, n, rp(:), ci(:)
    REAL(C_DOUBLE), INTENT(IN) :: val(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: tp(:), ti(:)
    REAL(C_DOUBLE), ALLOCATABLE, INTENT(OUT) :: tv(:)
    REAL(C_DOUBLE), OPTIONAL, INTENT(OUT) :: phase_times(4)
    REAL(C_DOUBLE) :: stamp,now
    INTEGER :: r, j, c, p
    INTEGER :: cursor(n)
    IF(PRESENT(phase_times)) stamp=CSRWallTime()
    ALLOCATE(tp(n+1), ti(rp(m+1)), tv(rp(m+1)))
    tp = 0
    IF(PRESENT(phase_times)) THEN
      now=CSRWallTime();phase_times(1)=now-stamp;stamp=now
    END IF
    DO j = 1, rp(m+1)
      c = ci(j)
      IF (c < 1 .OR. c > n) ERROR STOP 'TransposeCSR: invalid column'
      tp(c+1) = tp(c+1) + 1
    END DO
    IF(PRESENT(phase_times)) THEN
      now=CSRWallTime();phase_times(2)=now-stamp;stamp=now
    END IF
    DO c = 1, n
      tp(c+1) = tp(c+1) + tp(c)
    END DO
    cursor = tp(1:n)
    IF(PRESENT(phase_times)) THEN
      now=CSRWallTime();phase_times(3)=now-stamp;stamp=now
    END IF
    DO r = 1, m
      DO j = rp(r)+1, rp(r+1)
        c = ci(j)
        p = cursor(c)+1
        ti(p) = r
        tv(p) = val(j)
        cursor(c) = p
      END DO
    END DO
    IF(PRESENT(phase_times)) phase_times(4)=CSRWallTime()-stamp
  END SUBROUTINE

  ! Optional diagnostics use wall time without adding an MPI dependency.
  REAL(C_DOUBLE) FUNCTION CSRWallTime() RESULT(seconds)
    INTEGER(C_INT64_T) :: ticks,rate
    CALL SYSTEM_CLOCK(ticks,rate)
    seconds=REAL(ticks,C_DOUBLE)/REAL(rate,C_DOUBLE)
  END FUNCTION

  ! Bound each slab by its actual number of intermediate scalar products,
  ! rather than an independence/density model. A single k cannot be split
  ! further here; reject it if it exceeds the chosen workspace budget.
  SUBROUTINE PlanChunks(n, ci, brp, budget, ends, ok)
    INTEGER, INTENT(IN) :: n, ci(:), brp(:)
    INTEGER(C_INT64_T), INTENT(IN) :: budget
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ends(:)
    LOGICAL, INTENT(OUT) :: ok
    INTEGER(C_INT64_T) :: work(n), used
    INTEGER :: j, k, count, tmp(n+1)
    work = 0
    DO j = 1, SIZE(ci)
      k = ci(j)
      IF (k < 1 .OR. k > n) ERROR STOP 'PlanChunks: invalid column'
      work(k) = work(k) + INT(brp(k+1)-brp(k), C_INT64_T)
    END DO
    ok = ALL(work <= budget) .AND. budget > 0
    IF (.NOT. ok) RETURN
    count = 1
    tmp(1) = 0
    used = 0
    DO k = 1, n
      IF (used + work(k) > budget) THEN
        count = count+1
        tmp(count) = k-1
        used = 0
      END IF
      used = used+work(k)
    END DO
    count = count+1
    tmp(count) = n
    ends = tmp(1:count)
  END SUBROUTINE

  SUBROUTINE PackChunks(m, rp, ci, val, ends, crp, cci, cv, offsets)
    INTEGER, INTENT(IN) :: m, rp(:), ci(:), ends(:)
    REAL(C_DOUBLE), INTENT(IN) :: val(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: crp(:), cci(:), offsets(:)
    REAL(C_DOUBLE), ALLOCATABLE, INTENT(OUT) :: cv(:)
    INTEGER :: nc, ch, r, j, b, p
    INTEGER, ALLOCATABLE :: cursor(:)
    nc = SIZE(ends)-1
    ALLOCATE(crp(nc*(m+1)), cci(rp(m+1)), cv(rp(m+1)))
    ALLOCATE(offsets(nc+1), cursor(nc*m))
    crp = 0
    DO r = 1, m
      ch = 1
      DO j = rp(r)+1, rp(r+1)
        DO WHILE (ci(j) > ends(ch+1))
          ch = ch+1
        END DO
        b = (ch-1)*(m+1)
        crp(b+r+1) = crp(b+r+1)+1
      END DO
    END DO
    offsets(1) = 0
    DO ch = 1, nc
      b = (ch-1)*(m+1)
      DO r = 1, m
        crp(b+r+1) = crp(b+r+1)+crp(b+r)
      END DO
      offsets(ch+1) = offsets(ch)+crp(b+m+1)
    END DO
    cursor = 0
    DO r = 1, m
      ch = 1
      DO j = rp(r)+1, rp(r+1)
        DO WHILE (ci(j) > ends(ch+1))
          ch = ch+1
        END DO
        b = (ch-1)*(m+1)
        p = offsets(ch)+crp(b+r)+cursor((ch-1)*m+r)+1
        cci(p) = ci(j)-ends(ch)
        cv(p) = val(j)
        cursor((ch-1)*m+r) = cursor((ch-1)*m+r)+1
      END DO
    END DO
  END SUBROUTINE
END MODULE
