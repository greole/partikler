project(SPHIO LANGUAGES CXX)

add_library(SPHIO INTERFACE)
include(git_watcher.cmake)

target_include_directories(SPHIO INTERFACE ${PROJECT_SOURCE_DIR}/src)
