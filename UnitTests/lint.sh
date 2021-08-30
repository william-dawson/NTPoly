set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   # conda activate ntpoly-conda
   :
fi

# Python
flake8 UnitTests
flake8 Examples

# C/C++
clang-format -style=llvm Source/C/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.h -i
clang-format -style=llvm Source/CPlusPlus/*.cc -i
clang-format -style=llvm Examples/*/*.cc -i

# Fortran
for f in $(find . -type f -name "*.*90"); do
  emacs -batch $f --eval="(f90-mode)" \
                  --eval="(f90-upcase-keywords)" -f save-buffer 2>/dev/null
  emacs -batch $f --eval="(f90-mode)" \
                  --eval="(f90-indent-subprogram)" -f save-buffer 2>/dev/null
done

# 80 Column Limit
for f in $(find -L Source Targets UnitTests); do
	awk 'NF > 80 {print FILENAME ; print "Line " NR ; print ; stat = 1} \
	              END {exit stat}' FS= $f 2>/dev/null
done

# If git returns any changes, we know there is a problem.
git --no-pager diff --exit-code