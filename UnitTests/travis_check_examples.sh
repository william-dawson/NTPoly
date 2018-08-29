if [ -z ${TESTEXAMPLES+x} ]; then
  echo "Skipping examples"
else
  python testBuildInstructions.py ../Examples/ComplexMatrix/ run-fortran
  python testBuildInstructions.py ../Examples/ComplexMatrix/ run-c++
  python testBuildInstructions.py ../Examples/ComplexMatrix/ run-python
  python testBuildInstructions.py ../Examples/GraphTheory/ run-fortran
  python testBuildInstructions.py ../Examples/GraphTheory/ run-c++
  python testBuildInstructions.py ../Examples/GraphTheory/ run-python
  python testBuildInstructions.py ../Examples/HydrogenAtom/ run-fortran
  python testBuildInstructions.py ../Examples/HydrogenAtom/ run-c++
  python testBuildInstructions.py ../Examples/HydrogenAtom/ run-python
  python testBuildInstructions.py ../Examples/PremadeMatrix/ run-fortran
  python testBuildInstructions.py ../Examples/PremadeMatrix/ run-c++
  python testBuildInstructions.py ../Examples/PremadeMatrix/ run-python
fi
