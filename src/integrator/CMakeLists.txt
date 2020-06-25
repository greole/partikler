project(Integrator LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Integrator SHARED)

target_sources(Integrator
  PRIVATE
    src/PosIntegrator.cpp
  PUBLIC
    src/PosIntegrator.hpp
  )


add_dependencies(Integrator Core)
add_dependencies(Integrator Geometry)


target_link_libraries(Integrator
  PUBLIC
  Core
  Geometry
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES}
  )

target_include_directories(Integrator
  PUBLIC
  src)

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

#add_subdirectory(tests/SPHCore)
