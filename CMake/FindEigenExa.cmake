# Find the EigenExa Module
# Variables set:
# - EigenExa_FOUND - system found EigenExa.
# - EigenExa_LIBRARIES - the linker line for EigenExa.
# - EigenExa_INCLUDE_DIRS - the path to EigenExa.

# First we search for the libraries
find_library(EigenExa_LIBRARIES EigenExa)
find_path(EigenExa_INCLUDE_DIRS "eigen_libs_mod.mod")

# Now check if that worked
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EigenExa DEFAULT_MSG
  EigenExa_LIBRARIES EigenExa_INCLUDE_DIRS)
mark_as_advanced(EigenExa_LIBRARIES EigenExa_INCLUDE_DIRS)
