project(SPHKernel LANGUAGES CXX)

#find_package(GTest REQUIRED)
find_package(yaml-cpp REQUIRED)

add_library(SPHKernel SHARED
  src/STLWendland2D.cpp
  )

add_dependencies(SPHIO SPHIO)


target_link_libraries(SPHKernel SPHCore SPHIO OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(SPHKernel PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
