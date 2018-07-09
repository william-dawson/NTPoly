################################################################################
# Documentation
find_program (FORD ford)
find_package(Doxygen)

add_custom_target(doc)

# Fortran
if ( FORD )
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Ford.in
                 ${CMAKE_BINARY_DIR}/Documentation/Ford @ONLY)
  add_custom_target(fdoc "${FORD}" Ford)
  add_dependencies(doc fdoc)
else()
  MESSAGE("FORD not found! No fortran documentation will be generated.")
endif()

# C++
if (NOT FORTRAN_ONLY AND DOXYGEN_FOUND)
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Doxyfile.in
                 ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(cdoc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
  add_dependencies(doc cdoc)
elseif(DOXYGEN_FOUND)
  MESSAGE("Doxygen not found! No documentation will be generated.")
endif()
