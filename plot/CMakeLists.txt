add_library(
    plot
    plot.cpp
    plot.h
    matplotlibcpp.h
)

find_package(Python3 COMPONENTS Development NumPy)
target_include_directories(
    plot 
    PRIVATE 
    ${Python3_INCLUDE_DIRS} 
    ${Python3_NumPy_INCLUDE_DIRS} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../plot
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    plot
    src
    Python3::Python
    Python3::NumPy
)