project(SPHDatastructures LANGUAGES CXX)

add_library(SPHDatastructures INTERFACE)

target_include_directories(SPHDatastructures INTERFACE ${PROJECT_SOURCE_DIR}/src)
