!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE MatrixDRModule
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixSRModule, ONLY : Matrix_sr, ConstructMatrixSFromTripletList
  USE TripletListRModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList
  USE TripletRModule, ONLY : Triplet_r

#define DATATYPE REAL(NTREAL)
#define DMTYPE Matrix_dr
#define SMTYPE Matrix_sr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r

#include "includes/MatrixDImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixDRModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE MatrixDCModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixSCModule, ONLY : Matrix_sc, ConstructMatrixSFromTripletList
  USE TripletListCModule, ONLY : TripletList_c, ConstructTripletList, &
       & AppendToTripletList
  USE TripletCModule, ONLY : Triplet_c

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE Matrix_dc
#define SMTYPE Matrix_sc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c

#include "includes/MatrixDImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixDCModule
