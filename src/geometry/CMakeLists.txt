project(Geometry LANGUAGES CXX)

add_library(Geometry SHARED)

add_dependencies(Geometry IO Core)

target_sources(Geometry
  PRIVATE
    src/stl/STLLimitedDx.cpp
    src/Cubiod.cpp
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/stl/STLLimitedDx.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Cubiod.hpp"
  )

target_link_libraries(Geometry
  PUBLIC
    Core
  INTERFACE
    ${YAML_CPP_LIBRARIES}
  PRIVATE
    OpenMP::OpenMP_CXX
  )

target_include_directories(Geometry
  PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
