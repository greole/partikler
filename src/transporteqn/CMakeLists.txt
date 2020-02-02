project(TransportEqn LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(TransportEqn SHARED)

target_sources(TransportEqn
  PRIVATE
  "src/Momentum.cpp"
  "src/Conti.cpp"
  "src/Pressure.cpp"
  "src/Viscosity.cpp"
  "src/forces/Gravity.cpp"
  "src/forces/AdhesionAkinci.cpp"
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Momentum.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Conti.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Pressure.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Viscosity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/Gravity.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/forces/AdhesionAkinci.hpp"
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
)

#add_subdirectory(tests/SPHCore)
