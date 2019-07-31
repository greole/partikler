project(SPHCore LANGUAGES CXX)

find_package(GTest REQUIRED)
find_package(yaml-cpp REQUIRED)

add_library(SPHCore SHARED
  src/SPHModels.cpp
  src/Kernel.cpp
  src/ParticleGenerator.cpp
  src/ParticleNeighbours.cpp
  src/SortParticles.cpp
  src/Conti.cpp
  src/Pressure.cpp
  src/Momentum.cpp
  src/STLLimitedDx.cpp
  src/STLPosIntegrator.cpp
  src/Viscosity.cpp
  )

add_dependencies(SPHIO SPHIO)

target_link_libraries(SPHCore SPHIO OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARY})
target_include_directories(SPHCore PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests/SPHCore)
