add_library(
    src
    fmt.cpp
    fmt.h
    rrt.cpp
    rrt.h
    prm.cpp
    prm.h
    geometry.cpp
    geometry.h
    planning.cpp
    planning.h
)

target_include_directories(src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})