project(Geometry LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Geometry SHARED
  src/stl/STLLimitedDx.cpp
  src/stl/computes.cpp
  )

add_dependencies(Geometry Core)


target_link_libraries(Geometry Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
