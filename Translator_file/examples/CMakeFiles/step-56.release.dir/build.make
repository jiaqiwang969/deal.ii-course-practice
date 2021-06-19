# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /workspaces/dealii

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workspaces/dealii

# Include any dependencies generated for this target.
include examples/CMakeFiles/step-56.release.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/step-56.release.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/step-56.release.dir/flags.make

examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.o: examples/CMakeFiles/step-56.release.dir/flags.make
examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.o: examples/step-56/step-56.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.o"
	cd /workspaces/dealii/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/step-56.release.dir/step-56/step-56.cc.o -c /workspaces/dealii/examples/step-56/step-56.cc

examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/step-56.release.dir/step-56/step-56.cc.i"
	cd /workspaces/dealii/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/examples/step-56/step-56.cc > CMakeFiles/step-56.release.dir/step-56/step-56.cc.i

examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/step-56.release.dir/step-56/step-56.cc.s"
	cd /workspaces/dealii/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/examples/step-56/step-56.cc -o CMakeFiles/step-56.release.dir/step-56/step-56.cc.s

# Object files for target step-56.release
step__56_release_OBJECTS = \
"CMakeFiles/step-56.release.dir/step-56/step-56.cc.o"

# External object files for target step-56.release
step__56_release_EXTERNAL_OBJECTS =

bin/step-56.release: examples/CMakeFiles/step-56.release.dir/step-56/step-56.cc.o
bin/step-56.release: examples/CMakeFiles/step-56.release.dir/build.make
bin/step-56.release: lib/libdeal_II.so.10.0.0-pre
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libtbb.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libz.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_iostreams.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_serialization.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_system.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_thread.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_regex.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_chrono.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_date_time.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libboost_atomic.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libumfpack.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libcholmod.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libccolamd.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libcolamd.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libcamd.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libsuitesparseconfig.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libamd.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libmetis.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libarpack.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/liblapack.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libblas.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libassimp.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libgmsh.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libgsl.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libgslcblas.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libmuparser.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKBO.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKBool.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKBRep.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKernel.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKFeat.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKFillet.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKG2d.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKG3d.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKGeomAlgo.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKGeomBase.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKHLR.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKIGES.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKMath.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKMesh.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKOffset.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKPrim.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKShHealing.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKSTEP.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKSTEPAttr.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKSTEPBase.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKSTEP209.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKSTL.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKTopAlgo.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libTKXSBase.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libsundials_idas.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libsundials_arkode.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libsundials_kinsol.so
bin/step-56.release: /usr/lib/aarch64-linux-gnu/libsundials_nvecserial.so
bin/step-56.release: examples/CMakeFiles/step-56.release.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/step-56.release"
	cd /workspaces/dealii/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/step-56.release.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/step-56.release.dir/build: bin/step-56.release

.PHONY : examples/CMakeFiles/step-56.release.dir/build

examples/CMakeFiles/step-56.release.dir/clean:
	cd /workspaces/dealii/examples && $(CMAKE_COMMAND) -P CMakeFiles/step-56.release.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/step-56.release.dir/clean

examples/CMakeFiles/step-56.release.dir/depend:
	cd /workspaces/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/dealii /workspaces/dealii/examples /workspaces/dealii /workspaces/dealii/examples /workspaces/dealii/examples/CMakeFiles/step-56.release.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/step-56.release.dir/depend

