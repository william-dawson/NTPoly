set -e

# Python
flake8 UnitTests
flake8 Examples

# C/C++
clang-format -style=llvm Source/C/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.cc -i
clang-format -style=llvm Examples/*/*.cc -i

# If git returns any changes, we know there is a problem.
git --no-pager diff --exit-code

