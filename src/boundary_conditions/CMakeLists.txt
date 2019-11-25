project(BoundaryConditions LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(BoundaryConditions SHARED
  src/stl/GenerateBoundaryParticles.cpp
  )

add_dependencies(BoundaryConditions Core)


target_link_libraries(BoundaryConditions Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(BoundaryConditions PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
