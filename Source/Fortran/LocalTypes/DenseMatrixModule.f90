!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DenseMatrixModule
  USE DataTypesModule, ONLY : NTREAL
  USE SparseMatrixModule, ONLY : SparseMatrix_t, ConstructFromTripletList
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
       & AppendToTripletList
  USE TripletModule, ONLY : Triplet_t

#define DATATYPE REAL(NTREAL)
#define DMTYPE DenseMatrix_t
#define SMTYPE SparseMatrix_t
#define TTYPE Triplet_t
#define TLISTTYPE TripletList_t

#include "includes/DenseMatrixImplementation.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE DenseMatrixModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DenseMatrixCModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE SparseMatrixCModule, ONLY : SparseMatrix_c, ConstructFromTripletList
  USE TripletListCModule, ONLY : TripletList_c, ConstructTripletList, &
       & AppendToTripletList
  USE TripletCModule, ONLY : Triplet_c

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE DenseMatrix_c
#define SMTYPE SparseMatrix_c
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c

#include "includes/DenseMatrixImplementation.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE DenseMatrixCModule
