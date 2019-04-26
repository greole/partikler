project(SPHIO LANGUAGES CXX)

find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

include(git_watcher.cmake)

add_library(SPHIO SHARED src/RunTime)

add_dependencies(SPHIO SPHDatastructures)

target_link_libraries(SPHIO SPHDatastructures OpenMP::OpenMP_CXX) 

target_include_directories(SPHIO PUBLIC ${PROJECT_SOURCE_DIR}/src)
