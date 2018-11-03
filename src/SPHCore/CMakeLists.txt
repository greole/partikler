project(SPHCore LANGUAGES CXX)

find_package(GTest REQUIRED)
add_library(SPHCore INTERFACE)

target_include_directories(SPHCore INTERFACE ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests)
