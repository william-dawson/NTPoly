!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixSModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE MatrixMarketModule, ONLY : ParseMMHeader
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, real CSR matrix.
  TYPE, PUBLIC :: Matrix_lsr
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a local, complex CSR matrix.
  TYPE, PUBLIC :: Matrix_lsc
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Construct/Destruct
  PUBLIC :: DestructMatrix
  PUBLIC :: CopyMatrix
  !! Basic Accessors
  PUBLIC :: GetMatrixRows
  PUBLIC :: GetMatrixColumns
  PUBLIC :: ExtractMatrixRow
  PUBLIC :: ExtractMatrixColumn
  !! Routines for splitting and composing
  PUBLIC :: SplitMatrix
  PUBLIC :: SplitMatrixColumns
  PUBLIC :: ComposeMatrix
  PUBLIC :: ComposeMatrixColumns
  !! ETC
  PUBLIC :: TransposeMatrix
  PUBLIC :: PrintMatrix
  PUBLIC :: MatrixToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE Matrix_lsr
     MODULE PROCEDURE ConstructEmptyMatrix_lsr
     MODULE PROCEDURE ConstructMatrixFromFile_lsr
     MODULE PROCEDURE ConstructMatrixFromTripletList_lsr
  END INTERFACE
  INTERFACE Matrix_lsc
     MODULE PROCEDURE ConstructEmptyMatrix_lsc
     MODULE PROCEDURE ConstructMatrixFromFile_lsc
     MODULE PROCEDURE ConstructMatrixFromTripletList_lsc
  END INTERFACE
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_lsr
     MODULE PROCEDURE DestructMatrix_lsc
  END INTERFACE
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_lsr
     MODULE PROCEDURE CopyMatrix_lsc
  END INTERFACE
  INTERFACE GetMatrixRows
     MODULE PROCEDURE GetMatrixRows_lsr
     MODULE PROCEDURE GetMatrixRows_lsc
  END INTERFACE
  INTERFACE GetMatrixColumns
     MODULE PROCEDURE GetMatrixColumns_lsr
     MODULE PROCEDURE GetMatrixColumns_lsc
  END INTERFACE
  INTERFACE ExtractMatrixRow
     MODULE PROCEDURE ExtractMatrixRow_lsr
     MODULE PROCEDURE ExtractMatrixRow_lsc
  END INTERFACE
  INTERFACE ExtractMatrixColumn
     MODULE PROCEDURE ExtractMatrixColumn_lsr
     MODULE PROCEDURE ExtractMatrixColumn_lsc
  END INTERFACE
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_lsr
     MODULE PROCEDURE SplitMatrix_lsc
  END INTERFACE
  INTERFACE SplitMatrixColumns
     MODULE PROCEDURE SplitMatrixColumns_lsr
     MODULE PROCEDURE SplitMatrixColumns_lsc
  END INTERFACE
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_lsr
     MODULE PROCEDURE ComposeMatrix_lsc
  END INTERFACE
  INTERFACE ComposeMatrixColumns
     MODULE PROCEDURE ComposeMatrixColumns_lsr
     MODULE PROCEDURE ComposeMatrixColumns_lsc
  END INTERFACE
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_lsr
     MODULE PROCEDURE TransposeMatrix_lsc
  END INTERFACE
  INTERFACE PrintMatrix
     MODULE PROCEDURE PrintMatrix_lsr
     MODULE PROCEDURE PrintMatrix_lsc
  END INTERFACE
  INTERFACE MatrixToTripletList
     MODULE PROCEDURE MatrixToTripletList_lsr
     MODULE PROCEDURE MatrixToTripletList_lsc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define SMTYPE Matrix_lsr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r

#define ConstructEmptyMatrix ConstructEmptyMatrix_lsr
#define ConstructMatrixFromFile ConstructMatrixFromFile_lsr
#define ConstructZeroMatrix ConstructZeroMatrix_lsr
#define ConstructMatrixFromTripletList ConstructMatrixFromTripletList_lsr
#define DestructMatrix DestructMatrix_lsr
#define CopyMatrix CopyMatrix_lsr
#define GetMatrixRows GetMatrixRows_lsr
#define GetMatrixColumns GetMatrixColumns_lsr
#define ExtractMatrixRow ExtractMatrixRow_lsr
#define ExtractMatrixColumn ExtractMatrixColumn_lsr
#define SplitMatrix SplitMatrix_lsr
#define SplitMatrixColumns SplitMatrixColumns_lsr
#define ComposeMatrix ComposeMatrix_lsr
#define ComposeMatrixColumns ComposeMatrixColumns_lsr
#define TransposeMatrix TransposeMatrix_lsr
#define PrintMatrix PrintMatrix_lsr
#define MatrixToTripletList MatrixToTripletList_lsr

#include "includes/MatrixSImpl.F90"

#undef ConstructEmptyMatrix
#undef ConstructMatrixFromFile
#undef ConstructZeroMatrix
#undef ConstructMatrixFromTripletList
#undef DestructMatrix
#undef CopyMatrix
#undef GetMatrixRows
#undef GetMatrixColumns
#undef ExtractMatrixRow
#undef ExtractMatrixColumn
#undef SplitMatrix
#undef SplitMatrixColumns
#undef ComposeMatrix
#undef ComposeMatrixColumns
#undef TransposeMatrix
#undef PrintMatrix
#undef MatrixToTripletList

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define SMTYPE Matrix_lsc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c

#define ConstructEmptyMatrix ConstructEmptyMatrix_lsc
#define ConstructMatrixFromFile ConstructMatrixFromFile_lsc
#define ConstructZeroMatrix ConstructZeroMatrix_lsc
#define ConstructMatrixFromTripletList ConstructMatrixFromTripletList_lsc
#define DestructMatrix DestructMatrix_lsc
#define CopyMatrix CopyMatrix_lsc
#define GetMatrixRows GetMatrixRows_lsc
#define GetMatrixColumns GetMatrixColumns_lsc
#define ExtractMatrixRow ExtractMatrixRow_lsc
#define ExtractMatrixColumn ExtractMatrixColumn_lsc
#define SplitMatrix SplitMatrix_lsc
#define SplitMatrixColumns SplitMatrixColumns_lsc
#define ComposeMatrix ComposeMatrix_lsc
#define ComposeMatrixColumns ComposeMatrixColumns_lsc
#define TransposeMatrix TransposeMatrix_lsc
#define PrintMatrix PrintMatrix_lsc
#define MatrixToTripletList MatrixToTripletList_lsc

#define ISCOMPLEX 1
#include "includes/MatrixSImpl.F90"
#undef ISCOMPLEX

#undef ConstructEmptyMatrix
#undef ConstructMatrixFromFile
#undef ConstructZeroMatrix
#undef ConstructMatrixFromTripletList
#undef DestructMatrix
#undef CopyMatrix
#undef GetMatrixRows
#undef GetMatrixColumns
#undef ExtractMatrixRow
#undef ExtractMatrixColumn
#undef SplitMatrix
#undef SplitMatrixColumns
#undef ComposeMatrix
#undef ComposeMatrixColumns
#undef TransposeMatrix
#undef PrintMatrix
#undef MatrixToTripletList

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DATATYPE

END MODULE MatrixSModule
