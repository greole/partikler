project(Neighbours LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Neighbours SHARED)


target_sources(Neighbours
  PRIVATE
    src/ParticleNeighbours.cpp
    src/SortParticles.cpp
  PUBLIC
    src/ParticleNeighbours.hpp
    src/SortParticles.hpp
  )

add_dependencies(Neighbours Core)

target_link_libraries(Neighbours Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})

target_include_directories(Neighbours
  PUBLIC
  src
  )

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

#add_subdirectory(tests/SPHCore)
