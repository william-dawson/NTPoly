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

END MODULE TripletListCModule
