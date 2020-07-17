add_library(
    src
    fmt.cpp
    geometry.cpp
    fmt.h
    geometry.h
)

target_include_directories(src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})