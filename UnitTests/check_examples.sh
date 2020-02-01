if [ -z ${TESTEXAMPLES+x} ]; then
  echo "Skipping examples"
else
  python test_build.py ../Examples/ComplexMatrix/ run-fortran
  python test_build.py ../Examples/ComplexMatrix/ run-c++
  python test_build.py ../Examples/ComplexMatrix/ run-python
  python test_build.py ../Examples/GraphTheory/ run-fortran
  python test_build.py ../Examples/GraphTheory/ run-c++
  python test_build.py ../Examples/GraphTheory/ run-python
  python test_build.py ../Examples/HydrogenAtom/ run-fortran
  python test_build.py ../Examples/HydrogenAtom/ run-c++
  python test_build.py ../Examples/HydrogenAtom/ run-python
  python test_build.py ../Examples/MatrixMaps/ run-fortran
  python test_build.py ../Examples/MatrixMaps/ run-c++
  python test_build.py ../Examples/PremadeMatrix/ run-fortran
  python test_build.py ../Examples/PremadeMatrix/ run-c++
  python test_build.py ../Examples/PremadeMatrix/ run-python
fi
