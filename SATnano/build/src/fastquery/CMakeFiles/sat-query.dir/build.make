# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.6

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build

# Include any dependencies generated for this target.
include src/fastquery/CMakeFiles/sat-query.dir/depend.make

# Include the progress variables for this target.
include src/fastquery/CMakeFiles/sat-query.dir/progress.make

# Include the compile flags for this target's objects.
include src/fastquery/CMakeFiles/sat-query.dir/flags.make

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o: src/fastquery/CMakeFiles/sat-query.dir/flags.make
src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o: ../src/fastquery/main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o"
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/sat-query.dir/main.cpp.o -c /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/src/fastquery/main.cpp

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sat-query.dir/main.cpp.i"
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/src/fastquery/main.cpp > CMakeFiles/sat-query.dir/main.cpp.i

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sat-query.dir/main.cpp.s"
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/src/fastquery/main.cpp -o CMakeFiles/sat-query.dir/main.cpp.s

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.requires:
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.requires

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.provides: src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.requires
	$(MAKE) -f src/fastquery/CMakeFiles/sat-query.dir/build.make src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.provides.build
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.provides

src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.provides.build: src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.provides.build

# Object files for target sat-query
sat__query_OBJECTS = \
"CMakeFiles/sat-query.dir/main.cpp.o"

# External object files for target sat-query
sat__query_EXTERNAL_OBJECTS =

bin/sat-query: src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o
bin/sat-query: lib/libproc.a
bin/sat-query: lib/lib6mer.a
bin/sat-query: lib/libwavelib.a
bin/sat-query: lib/libutil.a
bin/sat-query: src/fastquery/CMakeFiles/sat-query.dir/build.make
bin/sat-query: src/fastquery/CMakeFiles/sat-query.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/sat-query"
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sat-query.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/fastquery/CMakeFiles/sat-query.dir/build: bin/sat-query
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/build

src/fastquery/CMakeFiles/sat-query.dir/requires: src/fastquery/CMakeFiles/sat-query.dir/main.cpp.o.requires
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/requires

src/fastquery/CMakeFiles/sat-query.dir/clean:
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery && $(CMAKE_COMMAND) -P CMakeFiles/sat-query.dir/cmake_clean.cmake
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/clean

src/fastquery/CMakeFiles/sat-query.dir/depend:
	cd /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/src/fastquery /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery /mnt/hgfs/data/projects/cwdtw-fastalign/data/cwSDTWnano/SATnano/build/src/fastquery/CMakeFiles/sat-query.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/fastquery/CMakeFiles/sat-query.dir/depend

