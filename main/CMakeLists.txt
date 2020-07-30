add_executable(
    main
    main.cpp
)

target_include_directories(
    main 
    PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../plot
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    main
    src
    plot
)