project(SPHDatastructures LANGUAGES CXX)

find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

add_library(SPHDatastructures SHARED src/SPHDatastructures)

target_include_directories(SPHDatastructures PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests)
