!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixSRModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL
  USE TripletListRModule, ONLY: TripletList_r, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletRModule, ONLY : Triplet_r

#define DATATYPE REAL(NTREAL)
#define SMTYPE Matrix_sr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r

#include "includes/MatrixSImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

END MODULE MatrixSRModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixSCModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX, MPINTREAL
  USE TripletListCModule, ONLY: TripletList_c, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletCModule, ONLY : Triplet_c

#define DATATYPE COMPLEX(NTCOMPLEX)
#define SMTYPE Matrix_sc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c

#include "includes/MatrixSImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

END MODULE MatrixSCModule
