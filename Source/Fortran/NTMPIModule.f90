!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module wraps the MPI include statement because on certain platforms
!! just writing "USE MPI" does not work.
MODULE NTMPIModule
#if USE_MPIF08
  USE MPI_F08
#elif USE_MPIH
  include "mpif.h"
#else
  USE MPI
#endif
END MODULE
