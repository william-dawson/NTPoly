################################################################################
# Determine python include path using python's distutils module
macro(get_py_include)
  execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c
      "from distutils.sysconfig import get_python_inc; print(get_python_inc())"
      OUTPUT_VARIABLE PYTHON_INCLUDE_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

################################################################################
# Determine Python Library Path
macro(get_py_lib)
  execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c
"
from distutils.sysconfig import get_python_lib
path=get_python_lib(standard_lib=True)+\"/../../Python\"
print(path)"
    OUTPUT_VARIABLE PYTHON_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()
