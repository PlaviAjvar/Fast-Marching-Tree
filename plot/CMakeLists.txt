add_library(
    plot
    plot.cpp
    plot.h
    matplotlibcpp.h
)

find_package(Python2 COMPONENTS Development NumPy)
target_include_directories(
    plot 
    PRIVATE 
    ${Python2_INCLUDE_DIRS} 
    ${Python2_NumPy_INCLUDE_DIRS} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../plot
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    plot
    src
    Python2::Python
    Python2::NumPy
)