!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletListRModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL
  USE TripletRModule, ONLY : Triplet_r, CompareTriplets

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL
#define TLISTTYPE TripletList_r
#define TTYPE Triplet_r

#include "includes/TripletListImpl.f90"

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletListRModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletListCModule
  USE DataTypesModule, ONLY: NTCOMPLEX, MPINTCOMPLEX
  USE TripletCModule, ONLY : Triplet_c, CompareTriplets

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX
#define TLISTTYPE TripletList_c
#define TTYPE Triplet_c

#include "includes/TripletListImpl.f90"

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletListCModule
