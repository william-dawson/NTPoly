##########################################################################
''' @package test_build.py
A script that tests the build instructions of the examples.
'''
from sys import argv
from subprocess import call
from os import environ, chdir
from os.path import join

def parse_command(fin, num_commands=1):
    '''
    Keep parsing a command which is wrapped in the github markdown style
    code blocks.
    '''
    ret_strings = []
    fin.readline()
    for i in range(0, num_commands):
        return_string = ""
        temp_string = fin.readline()
        while(temp_string != '\n' and temp_string != "\`\`\`"):
            return_string = return_string + temp_string
            temp_string = fin.readline()
        ret_strings.append(return_string)
    return ret_strings


##########################################################################
if __name__ == "__main__":
    check_directory = argv[1]
    check_command = argv[2]

    env_var = environ.copy()

    build_commands = []
    run_commands = []
    with open(join(check_directory, "ReadMe.md"), 'r') as fin:
        while True:
            line = fin.readline()
            if not line:
                break
            if check_command == "run-fortran":
                if line == "Fortran Build Instructions:\n":
                    build_commands = parse_command(fin)
                if line == "And then run with:\n":
                    run_commands = parse_command(fin)
            elif check_command == "run-c++":
                if line == "C++ Build Instructions:\n":
                    build_commands = parse_command(fin, 2)
                if line == "And then run with:\n":
                    run_commands = parse_command(fin)
                pass
            elif check_command == "run-python":
                env_var["PYTHONPATH"] = "../../Build/python"
                if line == "Run with python:\n":
                    run_commands = parse_command(fin)

    for i in range(0, len(build_commands)):
        build_commands[i] = [x for x in build_commands[i].split() if x != "\\"]
    for i in range(0, len(run_commands)):
        run_commands[i] = [x for x in run_commands[i].split() if x != "\\"]

    chdir(check_directory)
    for bc in build_commands:
        check = call(bc)
        if check != 0:
            print("Build Error")
            exit(-1)

    for rc in run_commands:
        rc = " ".join(rc)
        check = call(rc, shell=True, env=env_var)
        if check != 0:
            print("Run Error")
            exit(-1)
