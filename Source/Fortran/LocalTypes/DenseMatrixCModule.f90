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

END MODULE DenseMatrixCModule
