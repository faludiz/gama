set(TEST_BASE_DIR ${PROJECT_SOURCE_DIR}/tests/gama-g3)
set(GAMA_LOCAL ${CMAKE_BINARY_DIR}/gama-g3)
#set(INPUT_DIR ${TEST_BASE_DIR}/input)
#set(INPUT_FILES )

#set(RESULT_DIR ${CMAKE_BINARY_DIR}/tests/gama-g3/results/${PROJECT_VERSION})
#file(MAKE_DIRECTORY ${RESULT_DIR})
#file(MAKE_DIRECTORY ${RESULT_DIR}/gama-g3)

# ------------------------------------------------------------------------
#
# check transformation XYZ <--> BLH
#
add_executable(check_ellipsoid_xyz2blh
    src/check_ellipsoid_xyz2blh.cpp $<TARGET_OBJECTS:libgama> )
add_executable(check_ellipsoid_xyz2blh_list
    src/check_ellipsoid_xyz2blh_list.cpp $<TARGET_OBJECTS:libgama> )
add_executable(geng3test
    src/geng3test.cpp $<TARGET_OBJECTS:libgama>
    src/geng3test-md.h )

# target_compile_definitions(geng3test PRIVATE GenG3_DEBUG)

add_test(NAME gama-g3-ellipsoid-xyz2blh COMMAND check_ellipsoid_xyz2blh)
add_test(NAME gama-g3-ellipsoid-xyz2blh_list
    COMMAND check_ellipsoid_xyz2blh_list)
