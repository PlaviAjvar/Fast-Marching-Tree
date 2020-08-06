cmake_minimum_required(VERSION 3.16)
set(fmt Fast-Marching-Tree)

project(fmt C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

add_subdirectory(googletest ${CMAKE_CURRENT_SOURCE_DIR}/build/googletest)
add_subdirectory(sources ${CMAKE_CURRENT_SOURCE_DIR}/build/sources)
add_subdirectory(main ${CMAKE_CURRENT_SOURCE_DIR}/build/main)
add_subdirectory(plot ${CMAKE_CURRENT_SOURCE_DIR}/build/plot)
add_subdirectory(unit ${CMAKE_CURRENT_SOURCE_DIR}/build/unit)