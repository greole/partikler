project(SPHCore LANGUAGES CXX)

add_library(SPHCore INTERFACE)

target_include_directories(SPHCore INTERFACE ${PROJECT_SOURCE_DIR}/src)
