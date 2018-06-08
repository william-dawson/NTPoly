!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module That Contains a Description of a Data Type
MODULE TypeDescriptorModule
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(C)
    !> Real Value Data
    ENUMERATOR :: DESC_REAL = 1
    !> Complex Value Data
    ENUMERATOR :: DESC_COMPLEX = 2
    !> Matrix is stored locally
    ENUMERATOR :: DESC_LOCAL = 1
    !> Matrix is distributed
    ENUMERATOR :: DESC_DISTRIBUTED = 2
    !> Matrix is sparse.
    ENUMERATOR :: DESC_SPARSE = 1
    !> Matrix is dense.
    ENUMERATOR :: DESC_DENSE = 2
    !> Matrix is a general matrix
    ENUMERATOR :: DESC_GENERAL =1
    !> Matrix is symmetric (or Hermitian)
    ENUMERATOR :: DESC_SYMMETRIC =2
    !> Matrix is skew symmetric.
    ENUMERATOR :: DESC_SKEW_SYMMETRIC=3
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TypeDescriptorModule
