project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

 find_package(CGAL QUIET)
# include( ${CGAL_USE_FILE} )

set(test_sources SearchCubes.cpp)

add_executable(Partikler_Test_SearchCubes ${test_sources})

target_link_libraries(Partikler_Test_SearchCubes Core ${GTEST_LIBRARIES} gtest gtest_main)

gtest_discover_tests(Partikler_Test_SearchCubes)

add_test(NAME Partikler_Test_SearchCubes COMMAND Partikler_Test_SearchCubes)
