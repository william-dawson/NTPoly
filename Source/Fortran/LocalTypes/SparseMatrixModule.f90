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

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

END MODULE SparseMatrixModule

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

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

END MODULE SparseMatrixCModule
