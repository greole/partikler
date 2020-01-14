project(Geometry LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(Geometry SHARED
  # src/stl/computes.cpp
  # $<TARGET_OBJECTS:Core>

  )

add_dependencies(Geometry Core)

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
