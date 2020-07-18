add_executable(
    test
    test.cpp
    matplotlibcpp.h
)

find_package(Python2 COMPONENTS Development NumPy)
target_include_directories(
    test 
    PRIVATE 
    ${Python2_INCLUDE_DIRS} 
    ${Python2_NumPy_INCLUDE_DIRS} 
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    test
    src
    Python2::Python
    Python2::NumPy
)