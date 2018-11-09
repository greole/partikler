project( SPHDatastructuresTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL QUIET)
include( ${CGAL_USE_FILE} )

set(test_sources SPHDatastructures.cpp)

include_directories (BEFORE "${SPHDatastructures_SOURCE_DIR}")

add_executable(sphdatastructures ${test_sources})

target_link_libraries(sphdatastructures ${GTEST_LIBRARIES} pthread SPHDatastructures)

gtest_discover_tests(sphdatastructures)
