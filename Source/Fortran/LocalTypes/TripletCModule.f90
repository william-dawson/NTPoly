!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Complex.
MODULE TripletCModule
  USE DataTypesModule, ONLY: NTCOMPLEX, MPINTCOMPLEX

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPIDATATYPE MPINTCOMPLEX
#define INNERTYPE Triplet_c
#include "includes/TripletImplementation.f90"

END MODULE TripletCModule
