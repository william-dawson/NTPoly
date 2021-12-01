!! Allocate Eigen Exa
CALL InitializeEigenExa(A, exa)

!! Allocate Memory
ALLOCATE(AD(exa%local_rows, exa%local_cols))
AD = 0
ALLOCATE(VD(exa%local_rows, exa%local_cols))
VD = 0
ALLOCATE(WD(exa%mat_dim))
WD = 0

!! Convert to EigenExa
CALL StartTimer("NTToEigen")
#ifdef ISCOMPLEX
   CALL NTToEigen_c(A, AD, exa)
#else
   CALL NTToEigen_r(A, AD, exa)
#endif
CALL StopTimer("NTToEigen")

!! Calculate
CALL StartTimer("EigenExaCompute")
#ifdef ISCOMPLEX
    CALL Compute_c(AD, VD, WD, exa)
#else
    CALL Compute_r(AD, VD, WD, exa)
#endif
CALL StopTimer("EigenExaCompute")

!! Convert Back
CALL StartTimer("EigenToNT")

CALL ConstructEmptyMatrix(eigenvectors, A)
#ifdef ISCOMPLEX
    CALL EigenToNT_c(VD, eigenvectors, solver_parameters, exa)
#else
    CALL EigenToNT_r(VD, eigenvectors, solver_parameters, exa)
#endif

CALL ConstructEmptyMatrix(eigenvalues, A)
CALL ExtractEigenvalues(WD, eigenvalues, exa)

CALL StopTimer("EigenToNT")

!! Cleanup
#ifdef ISCOMPLEX
    CALL CleanUp_c(AD, VD, WD)
#else
    CALL CleanUp_r(AD, VD, WD)
#endif

IF (solver_parameters%be_verbose) THEN
 CALL PrintMatrixInformation(eigenvectors)
 CALL ExitSubLog
END IF
