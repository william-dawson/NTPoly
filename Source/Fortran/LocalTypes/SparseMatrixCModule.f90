!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE SparseMatrixCModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX, MPINTREAL
  USE TripletListCModule, ONLY: TripletList_c, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletCModule, ONLY : Triplet_c

#define DATATYPE COMPLEX(NTCOMPLEX)
#define SMTYPE SparseMatrix_c
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c
#include "includes/SparseMatrixImplementation.f90"

END MODULE SparseMatrixCModule
