flake8 UnitTests
flake8 Examples
clang-format -style=llvm Source/C/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.cc -i
clang-format -style=llvm Examples/*/*.cc -i
git --no-pager diff
