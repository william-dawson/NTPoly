!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE SparseMatrixModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL
  USE TripletListModule, ONLY: TripletList_t, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletModule, ONLY : Triplet_t

#define DATATYPE REAL(NTREAL)
#define SMTYPE SparseMatrix_t
#define TTYPE Triplet_t
#define TLISTTYPE TripletList_t
#include "includes/SparseMatrixImplementation.f90"

END MODULE SparseMatrixModule
