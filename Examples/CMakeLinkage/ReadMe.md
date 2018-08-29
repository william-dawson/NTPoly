# CMake Linkage

This example describes how to link against NTPoly using CMAKE. All of the
source files are the same as in the `Examples/PremadeMatrix` directory.
Look at the `CMakeLists.txt` file for all of the details. This file
assumes you have already built NTPoly and installed it:

> make install

though of course you can set the path to the config files yourself. To test
this example out, create a build directory, and then type:

> cmake ..

In this CMake file, I've written the openmp and mpi parts in a more vanilla
way, so this might not work for cross compiling. But the main way of
getting NTPoly to link will be the same.

> find_package(NTPoly REQUIRED)

> target_link_libraries(my_executable NTPoly::NTPoly)

And for C++:
> find_package(NTPolyWrapper REQUIRED)

> find_package(NTPolyCPP REQUIRED)

> target_link_libraries(my_executable NTPoly::NTPolyWrapper NTPoly::NTPolyCPP)
