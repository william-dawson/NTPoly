%module NTPolySwig
%include "typemaps.i"
%apply double& OUTPUT { double& chemical_potential_out };
%apply double& OUTPUT { double& energy_value_out };
%apply double *OUTPUT { double *max_power_eig };
%apply double *OUTPUT { double *min_ger_eig };
%apply double *OUTPUT { double *max_ger_eig };
%{
#include "SolverBase.h"
#include "ChebyshevSolvers.h"
#include "DensityMatrixSolvers.h"
#include "EigenBounds.h"
#include "ExponentialSolvers.h"
#include "GeometryOptimization.h"
#include "HermiteSolvers.h"
#include "InverseSolvers.h"
#include "LinearSolvers.h"
#include "LoadBalancer.h"
#include "MatrixMemoryPool.h"
#include "MinimizerSolvers.h"
#include "Permutation.h"
#include "PMatrixMemoryPool.h"
#include "PSMatrix.h"
#include "Polynomial.h"
#include "ProcessGrid.h"
#include "RootSolvers.h"
#include "SignSolvers.h"
#include "SMatrix.h"
#include "SolverParameters.h"
#include "SquareRootSolvers.h"
#include "TrigonometrySolvers.h"
#include "Triplet.h"
#include "TripletList.h"
#include <complex>
using namespace NTPoly;
%}

%include <complex.i>
%include "std_string.i"

%include "SolverBase.h"
%include "ChebyshevSolvers.h"
%include "DensityMatrixSolvers.h"
%include "EigenBounds.h"
%include "ExponentialSolvers.h"
%include "GeometryOptimization.h"
%include "HermiteSolvers.h"
%include "InverseSolvers.h"
%include "LinearSolvers.h"
%include "LoadBalancer.h"
%include "MatrixMemoryPool.h"
%include "MinimizerSolvers.h"
%include "Permutation.h"
%include "Polynomial.h"
%include "ProcessGrid.h"
%include "PMatrixMemoryPool.h"
%include "PSMatrix.h"
%include "RootSolvers.h"
%include "SignSolvers.h"
%include "SMatrix.h"
%include "SolverParameters.h"
%include "SquareRootSolvers.h"
%include "TrigonometrySolvers.h"
%include "Triplet.h"
%include "TripletList.h"
