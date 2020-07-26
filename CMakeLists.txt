cmake_minimum_required(VERSION 3.10)
set(fmt Fast-Marching-Tree)

project(fmt C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(sources ${CMAKE_CURRENT_SOURCE_DIR}/build/sources)
add_subdirectory(test ${CMAKE_CURRENT_SOURCE_DIR}/build/test)
add_subdirectory(plot ${CMAKE_CURRENT_SOURCE_DIR}/build/plot)