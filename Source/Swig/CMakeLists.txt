################################################################################
if (SWIG_FOUND)
  # Probe some information about python
  if(NOT DEFINED PYTHON_EXECUTABLE)
    find_package(PythonInterp REQUIRED)
  endif()
  include(ConfigPython.cmake)
  if(APPLE)
    get_py_lib()
    MESSAGE(STATUS "Using Python Library Path:" ${PYTHON_LIBRARIES})
  endif()
  get_py_include()

  MESSAGE(STATUS "Using Python Include Path:" ${PYTHON_INCLUDE_PATH})

  # Link everything
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Source/C)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/Source/CPlusPlus)

  # SWIG Setup
  INCLUDE(${SWIG_USE_FILE})
  SET(CMAKE_SWIG_FLAGS "")
  set(Swigsrc
    NTPolySwig.i
  )
  foreach(file ${Swigsrc})
    SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES CPLUSPLUS ON)
  endforeach(file)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_BINARY_DIR}/python)
  SWIG_ADD_MODULE(NTPolySwig python ${Swigsrc})
  SWIG_LINK_LIBRARIES(NTPolySwig NTPolyCPP NTPolyWrapper NTPoly
                      ${PYTHON_LIBRARIES} ${TOOLCHAIN_LIBS})
  set_target_properties(_NTPolySwig PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
    LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/python
  )

else()
  MESSAGE(WARNING "Swig not found! No python bindings will be generated.")
endif()
