################################################################################
from sys import argv

################################################################################
if __name__ == "__main__":
    check_directory = argv[1]
    check_command = argv[2]

    with open(check_command, 'r') fin:
        for line in fin:
            if check_command == "run-fortran":
                if line == "Fortran Build Instructions:":
                    pass
                if line == "And then run with:":
                    pass
                pass
            elif check_command == "run-c++":
                if line == "C++ Build Instructions:":
                    pass
                if line == "And then run with:":
                    pass
                pass
            elif check_command == "run-python":
                if line == "Python version:":
                    pass
                if line == "Python version:":
                    pass
                pass
