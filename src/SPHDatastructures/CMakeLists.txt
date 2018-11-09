project(SPHDatastructures LANGUAGES CXX)

find_package(GTest REQUIRED)
add_library(SPHDatastructures INTERFACE)

target_include_directories(SPHDatastructures INTERFACE ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests)
