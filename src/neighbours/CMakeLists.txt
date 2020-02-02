project(Neighbours LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Neighbours SHARED)


target_sources(Neighbours
  PRIVATE
    src/stl/STLParticleNeighbours.cpp
    src/ParticleNeighbours.cpp
    src/SortParticles.cpp
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/stl/STLParticleNeighbours.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleNeighbours.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/SortParticles.hpp"
  )

add_dependencies(Neighbours Core)

target_link_libraries(Neighbours Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(Neighbours PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
