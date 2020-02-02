project(Integrator LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Integrator SHARED)

target_sources(Integrator
  PRIVATE
    "src/stl/STLPosIntegrator.cpp"
    "src/PosIntegrator.cpp"
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/stl/STLPosIntegrator.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/PosIntegrator.hpp"
  )


add_dependencies(Integrator Core)
add_dependencies(Integrator Geometry)


target_link_libraries(Integrator
  Core
  Geometry
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES}
  )

target_include_directories(Integrator PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
