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
include plot/CMakeFiles/plot.dir/depend.make

# Include the progress variables for this target.
include plot/CMakeFiles/plot.dir/progress.make

# Include the compile flags for this target's objects.
include plot/CMakeFiles/plot.dir/flags.make

plot/CMakeFiles/plot.dir/plot.cpp.o: plot/CMakeFiles/plot.dir/flags.make
plot/CMakeFiles/plot.dir/plot.cpp.o: ../plot/plot.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tanja.ah/Desktop/Fast-Marching-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object plot/CMakeFiles/plot.dir/plot.cpp.o"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && /usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/plot.dir/plot.cpp.o -c /Users/tanja.ah/Desktop/Fast-Marching-Tree/plot/plot.cpp

plot/CMakeFiles/plot.dir/plot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/plot.dir/plot.cpp.i"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tanja.ah/Desktop/Fast-Marching-Tree/plot/plot.cpp > CMakeFiles/plot.dir/plot.cpp.i

plot/CMakeFiles/plot.dir/plot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/plot.dir/plot.cpp.s"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tanja.ah/Desktop/Fast-Marching-Tree/plot/plot.cpp -o CMakeFiles/plot.dir/plot.cpp.s

# Object files for target plot
plot_OBJECTS = \
"CMakeFiles/plot.dir/plot.cpp.o"

# External object files for target plot
plot_EXTERNAL_OBJECTS =

plot/libplot.a: plot/CMakeFiles/plot.dir/plot.cpp.o
plot/libplot.a: plot/CMakeFiles/plot.dir/build.make
plot/libplot.a: plot/CMakeFiles/plot.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tanja.ah/Desktop/Fast-Marching-Tree/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libplot.a"
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && $(CMAKE_COMMAND) -P CMakeFiles/plot.dir/cmake_clean_target.cmake
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/plot.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
plot/CMakeFiles/plot.dir/build: plot/libplot.a

.PHONY : plot/CMakeFiles/plot.dir/build

plot/CMakeFiles/plot.dir/clean:
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot && $(CMAKE_COMMAND) -P CMakeFiles/plot.dir/cmake_clean.cmake
.PHONY : plot/CMakeFiles/plot.dir/clean

plot/CMakeFiles/plot.dir/depend:
	cd /Users/tanja.ah/Desktop/Fast-Marching-Tree/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tanja.ah/Desktop/Fast-Marching-Tree /Users/tanja.ah/Desktop/Fast-Marching-Tree/plot /Users/tanja.ah/Desktop/Fast-Marching-Tree/build /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot /Users/tanja.ah/Desktop/Fast-Marching-Tree/build/plot/CMakeFiles/plot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : plot/CMakeFiles/plot.dir/depend

