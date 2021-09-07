"""
Generate the source code to be included in ELSI.

The steps are:
    1) Run CMake with the -S -save-temps options.
    2) Make the Fortran library.
    3) Copy the preprocessed files into an exportable directory.
    4) Place the custom Header switch files.
    5) Lint.
"""
from os import mkdir, getcwd, chdir, system
from os.path import join, basename
from glob import glob


cmake_line = '''cmake .. -DFORTRAN_ONLY=Yes \
-DCMAKE_Fortran_FLAGS_RELEASE="-cpp -S -save-temps"'''

switch1 = '''MODULE NTMPIModule

  USE MPI

  IMPLICIT NONE

  PUBLIC

END MODULE NTMPIModule'''

switch2 = '''MODULE NTMPIModule

  IMPLICIT NONE

  INCLUDE "mpif.h"

  PUBLIC

END MODULE NTMPIModule'''

if __name__ == "__main__":
    # Create the build directory
    try:
        mkdir(join("ElsiBuild"))
    except IOError:
        pass

    # Create the output directory
    try:
        mkdir(join("ElsiBuild", "elsi_output"))
    except IOError:
        pass

    # Run cmake.
    start_dir = getcwd()
    try:
        chdir("ElsiBuild")
        system(cmake_line)
        system("make NTPoly")
    finally:
        chdir(start_dir)

    # Copy the files.
    file_list = glob(join("ElsiBuild", "Source", "Fortran", "CMakeFiles",
                          "NTPoly.dir", "*.f90"))
    for f in file_list:
        new_name = basename(f.replace(".F90", ""))
        with open(f) as ifile:
            source = "".join([x for x in ifile if x[0] != "#"])
        with open(join("ElsiBuild", "elsi_output", new_name), "w") as ofile:
            ofile.write(source)

    # Create the special MPI switches
    with open(join("ElsiBuild", "elsi_output",
                   "NTMPIModule.f90"), "w") as ofile:
        ofile.write(switch1)
    with open(join("ElsiBuild", "elsi_output",
                   "NTMPIFH.f90"), "w") as ofile:
        ofile.write(switch2)

    # Might as well lint the result as well.
    finished_list = glob(join("ElsiBuild", "elsi_output", "*.f90"))
    args = '''--eval="(setq make-backup-files nil)" ''' + \
           '''--eval="(f90-mode)" --eval="(f90-indent-subprogram)" ''' + \
           '''-f save-buffer'''
    for f in finished_list:
        system("emacs -batch " + f + " " + args + " 2>/dev/null")
