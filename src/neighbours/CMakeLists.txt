project(Neighbours LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Neighbours SHARED
  src/ParticleNeighbours.cpp
  src/SortParticles.cpp
  )

add_dependencies(Neighbours Core)


target_link_libraries(Neighbours Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(Neighbours PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
