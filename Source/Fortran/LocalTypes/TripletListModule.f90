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

END MODULE TripletListModule
