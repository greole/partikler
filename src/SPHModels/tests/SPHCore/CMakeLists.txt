project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL QUIET)
include( ${CGAL_USE_FILE} )

set(test_sources sph_core.cpp)

include_directories (BEFORE "${SPHCore_SOURCE_DIR}")

add_executable(sph_core ${test_sources})

target_link_libraries(sph_core ${GTEST_LIBRARIES} pthread SPHDatastructures SPHCore SPHIO)

gtest_discover_tests(sph_core)
