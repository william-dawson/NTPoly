%module NTPolySwig
%include "typemaps.i"
%apply double & INOUT { double & chemical_potential_out };
%{
#include "SolverBase.h"
#include "ChebyshevSolvers.h"
#include "DensityMatrixSolvers.h"
#include "DistributedBlockedSparseMatrix.h"
#include "DistributedMatrixMemoryPool.h"
#include "ExponentialSolvers.h"
#include "FixedSolversParameters.h"
#include "InverseSolvers.h"
#include "IterativeSolversParameters.h"
#include "LinearSolvers.h"
#include "LoadBalancer.h"
#include "MatrixMemoryPool.h"
#include "MinimizerSolvers.h"
#include "Permutation.h"
#include "Polynomial.h"
#include "ProcessGrid.h"
#include "RootSolvers.h"
#include "SignSolvers.h"
#include "SparseMatrix.h"
#include "SquareRootSolvers.h"
#include "TrigonometrySolvers.h"
#include "Triplet.h"
#include "TripletList.h"
using namespace NTPoly;
%}
%include "std_string.i"

%include "SolverBase.h"
%include "ChebyshevSolvers.h"
%include "DensityMatrixSolvers.h"
%include "DistributedBlockedSparseMatrix.h"
%include "DistributedMatrixMemoryPool.h"
%include "ExponentialSolvers.h"
%include "FixedSolversParameters.h"
%include "InverseSolvers.h"
%include "IterativeSolversParameters.h"
%include "LinearSolvers.h"
%include "LoadBalancer.h"
%include "MatrixMemoryPool.h"
%include "MinimizerSolvers.h"
%include "Permutation.h"
%include "Polynomial.h"
%include "ProcessGrid.h"
%include "RootSolvers.h"
%include "SignSolvers.h"
%include "SparseMatrix.h"
%include "SquareRootSolvers.h"
%include "TrigonometrySolvers.h"
%include "Triplet.h"
%include "TripletList.h"