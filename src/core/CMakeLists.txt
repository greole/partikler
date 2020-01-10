project(Core LANGUAGES CXX)

find_package(yaml-cpp REQUIRED)

# TODO move SearchCubes to neighbours?

add_library(Core SHARED)

target_sources(Core
  PRIVATE
    src/Vec3.cpp
    src/Object.cpp
    src/Field.cpp
    src/Models.cpp
    src/Helper.cpp
    src/SearchCubes.cpp
    src/ObjectRegistry.cpp
    src/cgal/CGALHelper.cpp
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Vec3.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Object.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Field.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Models.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Helper.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/SearchCubes.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ObjectRegistry.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cgal/CGALHelper.hpp"
  )


target_link_libraries(Core
  PUBLIC
    ${YAML_CPP_LIBRARIES}
    ${CGAL_LIBRARIES}
    OpenMP::OpenMP_CXX
  )

target_include_directories(Core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src
  )

add_subdirectory(tests)
