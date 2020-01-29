project(BoundaryConditions LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(BoundaryConditions SHARED)

target_sources(BoundaryConditions
  PRIVATE
    src/stl/GenerateBoundaryParticles.cpp
    src/FixedValue.cpp
    src/CreateFields.cpp
  PUBLIC
    src/stl/GenerateBoundaryParticles.hpp
    src/FixedValue.hpp
    src/CreateFields.hpp
  )

add_dependencies(BoundaryConditions
  Core
  Geometry
  )


target_link_libraries(BoundaryConditions
  PUBLIC
  Core
  Geometry
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES})

target_include_directories(
  BoundaryConditions
  PUBLIC ${PROJECT_SOURCE_DIR}/src
  )

#add_subdirectory(tests/SPHCore)
