################################################################################
# Probe some information about python
if(NOT DEFINED PYTHON_EXECUTABLE)
  find_package(PythonInterp REQUIRED)
endif()

include(ConfigPython.cmake)
get_py_include()
message(STATUS "Using Python Include Path:" ${PYTHON_INCLUDE_PATH})

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Source/C)
include_directories(${CMAKE_SOURCE_DIR}/Source/CPlusPlus)

include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "")

set(Swigsrc
  NTPolySwig.i
)

foreach(file ${Swigsrc})
  SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES CPLUSPLUS ON)
endforeach(file)
set(CMAKE_SWIG_OUTDIR ${CMAKE_BINARY_DIR}/python)
if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
   swig_add_module(NTPolySwig python ${Swigsrc})
else()
   swig_add_library(NTPolySwig LANGUAGE python SOURCES ${Swigsrc})
endif()
if (APPLE)
   set(APPLE_SUPRESS "-flat_namespace -undefined suppress")
endif()
swig_link_libraries(NTPolySwig NTPolyCPP NTPolyWrapper NTPoly 
                    ${EigenExa_LIBRARIES}
                    ${TOOLCHAIN_LIBS} ${APPLE_SUPRESS})
set_target_properties(_NTPolySwig PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
  LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/python
)
