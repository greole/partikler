project( PartiklerTests LANGUAGES CXX)

include_directories(${GTEST_SOURCE_DIR}/include)
include_directories(${GTEST_INCLUDE_DIRS})

find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)
include( ${CGAL_USE_FILE} )

set(CMAKE_INSTALL_RPATH "/usr/local/lib64")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(test_sources search_cubes_tests.cpp)

# include_directories (BEFORE "${SPHCore_SOURCE_DIR}")

add_executable(search_cubes_tests ${test_sources})

target_link_libraries(search_cubes_tests ${GTEST_LIBRARIES} CGAL pthread SPHDatastructures OpenMP::OpenMP_CXX)
target_compile_options(search_cubes_tests PRIVATE ${OpenMP_CXX_FLAGS})

gtest_discover_tests(search_cubes_tests)

add_executable(search_cubes_benchmark "search_cubes_benchmark.cpp")
target_link_libraries(search_cubes_benchmark benchmark pthread CGAL SPHDatastructures OpenMP::OpenMP_CXX)
