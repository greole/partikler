project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL QUIET)
include( ${CGAL_USE_FILE} )

set(test_sources search_cubes.cpp)

# include_directories (BEFORE "${SPHCore_SOURCE_DIR}")

add_executable(search_cubes ${test_sources})

target_link_libraries(search_cubes ${GTEST_LIBRARIES} pthread SPHDatastructures)

gtest_discover_tests(search_cubes)
