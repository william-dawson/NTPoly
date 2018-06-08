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

END MODULE DenseMatrixModule
