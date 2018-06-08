!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletRModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL
#define INNERTYPE Triplet_r

#include "includes/TripletImpl.f90"

#undef INNERTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletRModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Complex.
MODULE TripletCModule
  USE DataTypesModule, ONLY: NTCOMPLEX, MPINTCOMPLEX

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX
#define INNERTYPE Triplet_c

#include "includes/TripletImpl.f90"

#undef INNERTYPE
#undef MPIDATATYPE
#undef DATATYPE

END MODULE TripletCModule
