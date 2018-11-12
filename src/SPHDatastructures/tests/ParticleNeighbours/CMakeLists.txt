project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL QUIET)
include( ${CGAL_USE_FILE} )

set(test_sources particle_neighbours.cpp)

# include_directories (BEFORE "${SPHCore_SOURCE_DIR}")

add_executable(particle_neighbours ${test_sources})

target_link_libraries(particle_neighbours ${GTEST_LIBRARIES} pthread SPHDatastructures)

gtest_discover_tests(particle_neighbours)
