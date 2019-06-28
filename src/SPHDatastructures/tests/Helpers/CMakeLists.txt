project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL QUIET)
include( ${CGAL_USE_FILE} )

set(test_sources Helper_tests.cpp)

add_executable(helper_test ${test_sources})

target_link_libraries(helper_test ${GTEST_LIBRARIES} pthread SPHDatastructures)

gtest_discover_tests(helper_test)
