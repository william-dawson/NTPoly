"""
Some helper functions for the unit tests.
"""
from os import environ
from os.path import join
# Path to the scratch directory.
scratch_dir = environ['SCRATCHDIR']
# Filename to use for results.
result_file = join(scratch_dir, "result.mtx")
# Logfile for testing the logger
log_file = join(scratch_dir, "log.yaml")
# Threshold value for comparing floating point values.
THRESHOLD = 1e-4
# Threshold used for checking extrapolazation.
EXTRAPTHRESHOLD = 1e-1
