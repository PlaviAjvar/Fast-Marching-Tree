# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tanja.ah/Desktop/Fast-Marching-Tree

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tanja.ah/Desktop/Fast-Marching-Tree/build

# Include any dependencies generated for this target.
include sources/CMakeFiles/src.dir/depend.make

# Include the progress variables for this target.
include sources/CMakeFiles/src.dir/progress.make

# Include the compile flags for this target's objects.
include sources/CMakeFiles/src.dir/flags.make

sources/CMakeFiles/src.dir/fmt.cpp.o: sources/CMakeFiles/src.dir/flags.make
sources/CMakeFiles/src.dir/fmt.cpp.o: ../sources/fmt.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tanja.ah/Desktop/Fast-Marching-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object sources/CMakeFiles/src.dir/fmt.cpp.o"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/src.dir/fmt.cpp.o -c /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/fmt.cpp

sources/CMakeFiles/src.dir/fmt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/src.dir/fmt.cpp.i"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/fmt.cpp > CMakeFiles/src.dir/fmt.cpp.i

sources/CMakeFiles/src.dir/fmt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/src.dir/fmt.cpp.s"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/fmt.cpp -o CMakeFiles/src.dir/fmt.cpp.s

sources/CMakeFiles/src.dir/geometry.cpp.o: sources/CMakeFiles/src.dir/flags.make
sources/CMakeFiles/src.dir/geometry.cpp.o: ../sources/geometry.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tanja.ah/Desktop/Fast-Marching-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object sources/CMakeFiles/src.dir/geometry.cpp.o"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/src.dir/geometry.cpp.o -c /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/geometry.cpp

sources/CMakeFiles/src.dir/geometry.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/src.dir/geometry.cpp.i"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/geometry.cpp > CMakeFiles/src.dir/geometry.cpp.i

sources/CMakeFiles/src.dir/geometry.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/src.dir/geometry.cpp.s"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources/geometry.cpp -o CMakeFiles/src.dir/geometry.cpp.s

# Object files for target src
src_OBJECTS = \
"CMakeFiles/src.dir/fmt.cpp.o" \
"CMakeFiles/src.dir/geometry.cpp.o"

# External object files for target src
src_EXTERNAL_OBJECTS =

sources/libsrc.a: sources/CMakeFiles/src.dir/fmt.cpp.o
sources/libsrc.a: sources/CMakeFiles/src.dir/geometry.cpp.o
sources/libsrc.a: sources/CMakeFiles/src.dir/build.make
sources/libsrc.a: sources/CMakeFiles/src.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tanja.ah/Desktop/Fast-Marching-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libsrc.a"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && $(CMAKE_COMMAND) -P CMakeFiles/src.dir/cmake_clean_target.cmake
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/src.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
sources/CMakeFiles/src.dir/build: sources/libsrc.a

.PHONY : sources/CMakeFiles/src.dir/build

sources/CMakeFiles/src.dir/clean:
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources && $(CMAKE_COMMAND) -P CMakeFiles/src.dir/cmake_clean.cmake
.PHONY : sources/CMakeFiles/src.dir/clean

sources/CMakeFiles/src.dir/depend:
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tanja.ah/Desktop/Fast-Marching-Tree /Users/tanja.ah/Desktop/Fast-Marching-Tree/sources /Users/tanja.ah/Desktop/Fast-Marching-Tree/build /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/sources/CMakeFiles/src.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sources/CMakeFiles/src.dir/depend

