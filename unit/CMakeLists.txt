add_executable(
    unit
    unit.cpp
)

target_include_directories(
    unit
    PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../sources
    ${CMAKE_CURRENT_SOURCE_DIR}/../googletest
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    unit
    src
    gtest_main
)

add_test(
    NAME equals
    COMMAND equals
)

add_test(
    NAME normalize
    COMMAND normalize
)

add_test(
    NAME dotproduct
    COMMAND dotproduct
)

add_test(
    NAME crossproduct
    COMMAND crossproduct
)

add_test(
    NAME weighedEuclidean
    COMMAND weighedEuclidean
)

add_test(
    NAME online
    COMMAND online
)

add_test(
    NAME intersect2D
    COMMAND intersect2D
)

add_test(
    NAME linedistance
    COMMAND linedistance
)

add_test(
    NAME signedscaler
    COMMAND signedscaler
)

add_test(
    NAME intersectbox,
    COMMAND intersectbox
)

add_test(
    NAME intersect3D
    COMMAND intersect3D
)

add_test(
    NAME collides2D
    COMMAND collides2D
)

add_test(
    NAME collides3D
    COMMAND collides3D
)