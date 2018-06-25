%module NTPolySwig
%include "typemaps.i"
%apply double& OUTPUT { double& chemical_potential_out };
%apply double *OUTPUT { double *max_power_eig };
%apply double *OUTPUT { double *min_ger_eig };
%apply double *OUTPUT { double *max_ger_eig };
%{
#include "MatrixMemoryPool.h"
#include "ProcessGrid.h"
#include "SparseMatrix.h"
#include "TrigonometrySolvers.h"
#include "Triplet.h"
#include "TripletList.h"
#include <complex.h>
using namespace NTPoly;
%}

%include <complex.i>
%include "std_string.i"

%include "MatrixMemoryPool.h"
%include "ProcessGrid.h"
%include "SparseMatrix.h"
%include "Triplet.h"
%include "TripletList.h"
