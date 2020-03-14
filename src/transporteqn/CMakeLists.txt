project(TransportEqn LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(TransportEqn SHARED)

target_sources(TransportEqn
  PRIVATE
  "src/Momentum.cpp"
  "src/continuity/Conti.cpp"
  "src/viscosity/Szewc.cpp"
  "src/pressure/Bonet.cpp"
  "src/pressure/Solenthaler.cpp"
  "src/forces/ForcesList.cpp"
  "src/forces/Gravity.cpp"
  "src/forces/surface_tension/Akinci.cpp"
  "src/forces/surface_tension/AkinciCurvature.cpp"
  PUBLIC
  "src/Momentum.hpp"
  "src/continuity/Conti.hpp"
  "src/viscosity/Szewc.hpp"
  "src/pressure/Bonet.hpp"
  "src/pressure/Solenthaler.hpp"
  "src/forces/ForcesList.hpp"
  "src/forces/Gravity.hpp"
  "src/forces/surface_tension/Akinci.hpp"
  "src/forces/surface_tension/AkinciCurvature.hpp"
  )

add_dependencies(TransportEqn Core)

target_link_libraries(TransportEqn
  Core
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES}
)

target_include_directories(TransportEqn
  PUBLIC
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/src/forces
  ${PROJECT_SOURCE_DIR}/src/forces/surface_tension
  ${PROJECT_SOURCE_DIR}/src/pressure
  ${PROJECT_SOURCE_DIR}/src/viscosity
  ${PROJECT_SOURCE_DIR}/src/continuity
)

#add_subdirectory(tests/SPHCore)
