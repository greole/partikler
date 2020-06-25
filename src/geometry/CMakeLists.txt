project(Geometry LANGUAGES CXX)

add_library(Geometry SHARED)

add_dependencies(Geometry IO Core)

target_sources(Geometry
  PRIVATE
    src/ParticleGeneratorBase.cpp
    src/Cubiod.cpp
    src/SurfaceNormal.cpp
    src/SelectMaterial.cpp
  PUBLIC
    src/ParticleGeneratorBase.hpp
    src/Cubiod.hpp
    src/SurfaceNormal.hpp
    src/SelectMaterial.hpp
  )

target_link_libraries(Geometry
  PUBLIC
    Core
    IO
  INTERFACE
    ${YAML_CPP_LIBRARIES}
  PRIVATE
    OpenMP::OpenMP_CXX
  )

target_include_directories(Geometry
  PUBLIC
    src
  )

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

#add_subdirectory(tests/SPHCore)
