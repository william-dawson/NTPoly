!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL
  USE TripletModule, ONLY : Triplet_t, CompareTriplets

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL
#define TLISTTYPE TripletList_t
#define TTYPE Triplet_t

#include "includes/TripletListImplementation.f90"

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletListModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletListCModule
  USE DataTypesModule, ONLY: NTCOMPLEX, MPINTCOMPLEX
  USE TripletCModule, ONLY : Triplet_c, CompareTriplets

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX
#define TLISTTYPE TripletList_c
#define TTYPE Triplet_c

#include "includes/TripletListImplementation.f90"

#undef TTYPE
#undef TLISTTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletListCModule
