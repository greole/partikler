project(Integrator LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Integrator SHARED
  src/stl/STLPosIntegrator.cpp
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
