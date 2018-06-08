!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Triplets of Integer, Integer, Double.
MODULE TripletModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL

#define DATATYPE REAL(NTREAL)
#define MPIDATATYPE MPINTREAL
#define INNERTYPE Triplet_t
#include "includes/TripletImplementation.f90"

END MODULE TripletModule
