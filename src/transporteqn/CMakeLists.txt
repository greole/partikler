project(SPHTransportEqn LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(SPHTransportEqn SHARED
  src/Momentum.cpp
  src/Conti.cpp
  src/Pressure.cpp
  src/Viscosity.cpp
  )

add_dependencies(SPHIO SPHIO)


target_link_libraries(SPHTransportEqn SPHCore SPHIO OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(SPHTransportEqn PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
