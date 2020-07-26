add_executable(
    test
    test.cpp
)

target_include_directories(
    test 
    PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../plot
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    test
    src
    plot
)