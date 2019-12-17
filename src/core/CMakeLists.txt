project(Core LANGUAGES CXX)

find_package(yaml-cpp REQUIRED)

# TODO move SearchCubes to neighbours?

add_library(Core SHARED
  src/Vec3.cpp
  src/FileIO.cpp
  src/Object.cpp
  src/Field.cpp
  src/Models.cpp
  src/SearchCubes.cpp
  src/Datastructures.cpp
  src/ObjectRegistry.cpp
  src/cgal/CGALHelper.cpp
  )


target_link_libraries(Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES} ${libCGAL})
target_include_directories(Core PUBLIC ${PROJECT_SOURCE_DIR}/src ${cgal_SOURCE_DIR})

add_subdirectory(tests)
