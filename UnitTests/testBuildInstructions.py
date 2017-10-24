################################################################################
from sys import argv
import subprocess
import os

################################################################################
def parse_command(fin):
    return_string = ""
    temp_string = fin.readline()
    while(temp_string != '\n'):
        return_string = return_string + temp_string
        temp_string = fin.readline()
        print(":::",temp_string)
    return return_string

################################################################################
if __name__ == "__main__":
    check_directory = argv[1]
    check_command = argv[2]

    build_command = ""
    run_command = ""
    with open(check_directory+"/ReadMe.md", 'r') as fin:
        while True:
            line = fin.readline()
            if not line: break
            if check_command == "run-fortran":
                if line == "Fortran Build Instructions:\n":
                    build_command = parse_command(fin)
                if line == "And then run with:\n":
                    run_command = parse_command(fin)
            elif check_command == "run-c++":
                if line == "C++ Build Instructions:\n":
                    build_command = parse_command(fin)
                if line == "And then run with:\n":
                    run_command = parse_command(fin)
                pass
            elif check_command == "run-python":
                if line == "Setup python environment:\n":
                    build_command = parse_command(fin)
                if line == "Run with python:\n":
                    run_command = parse_command(fin)

    build_command = [x for x in build_command.split() if x != "\\" ]
    run_command = [x for x in run_command.split() if x != "\\" ]

    os.chdir(check_directory)
    if check_command != "run-python":
        subprocess.run(build_command, check=True)
        subprocess.run(run_command, check=True)
    else:
        print(run_command)
        env_var = os.environ.copy()
        env_var["PYTHONPATH"] = "../../Build/python"
        run_command = " ".join(run_command)
        subprocess.run(run_command, check=True, shell=True, env=env_var)
