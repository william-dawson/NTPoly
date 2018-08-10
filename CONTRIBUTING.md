# Contributing Guide

In this page, we would like to provide information for developers willing
to help contribute to this project.

## Developer Work Flow

Clone the repository and start a branch for whatever work you are doing. If
your project is well formulated, you can also raise it as an enhancement in
the issues tab, and we can assign the work to you.

## Coding Style

All public subroutines should be documented using Doxygen. Private routines can
also optionally be documented this way.

All text that is written to the console should be written through the logging
module (this makes it easier to parse the output).

The source code is beautified using the atom beautify package. If you don't
want to install a beautifier that is ok, just mention that you need someone
to run a beautify script in your pull request. The C++ beautifier is set to
clang-format. 

If your code makes use of something from the literature, put a reference to it
in the Citations.bib file.

For floating point values, please use the types in the DataTypesModule. This
helps maintain C compatibility and makes it easier to change the precision.
