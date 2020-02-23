project(TransportEqn LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(TransportEqn SHARED)

target_sources(TransportEqn
  PRIVATE
  "src/Momentum.cpp"
  "src/continuity/Conti.cpp"
  "src/pressure/Bonet.cpp"
  "src/Viscosity.cpp"
  "src/forces/Gravity.cpp"
  "src/forces/surface_tension/Akinci.cpp"
  "src/forces/surface_tension/AkinciCurvature.cpp"
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Momentum.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/continuity/Conti.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pressure/Bonet.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Viscosity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/Gravity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/surface_tension/Akinci.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/surface_tension/AkinciCurvature.hpp"
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
  ${PROJECT_SOURCE_DIR}/src/continuity
)

#add_subdirectory(tests/SPHCore)
