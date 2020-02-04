project(TransportEqn LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(TransportEqn SHARED)

target_sources(TransportEqn
  PRIVATE
  "src/Momentum.cpp"
  "src/Conti.cpp"
  "src/pressure/PressureBonet.cpp"
  "src/Viscosity.cpp"
  "src/forces/Gravity.cpp"
  "src/forces/surface_tension/Akinci.cpp"
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Momentum.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Conti.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/pressure/PressureBonet.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Viscosity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/Gravity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/surface_tension/Akinci.hpp"
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
)

#add_subdirectory(tests/SPHCore)
