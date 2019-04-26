project(SPHCore LANGUAGES CXX)

find_package(GTest REQUIRED)

add_library(SPHCore SHARED src/SPHCore src/computes)

add_dependencies(SPHIO SPHIO)

target_link_libraries(SPHCore SPHIO OpenMP::OpenMP_CXX) 
target_include_directories(SPHCore PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests/SPHCore)
