################################################################################
# Determine python include path using python's distutils module
macro(get_py_include)
  execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c
      "from distutils.sysconfig import get_python_inc; print(get_python_inc())"
      OUTPUT_VARIABLE PYTHON_INCLUDE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()
