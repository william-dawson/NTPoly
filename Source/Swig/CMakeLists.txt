################################################################################
if (SWIG_FOUND)
  # Probe some information about python
  if(NOT DEFINED PYTHON_EXECUTABLE)
    find_package(PythonInterp REQUIRED)
  endif()

  include(ConfigPython.cmake)
  get_py_lib()
  message(STATUS "Using Python Library Path:" ${PYTHON_LIBRARIES})
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
