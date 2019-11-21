project(Core LANGUAGES CXX)

#find_package(GTest REQUIRED)
find_package(yaml-cpp REQUIRED)

# TODO move SearchCubes to neighbours?

add_library(Core SHARED
  src/FileIO.cpp
  src/Object.cpp
  src/Models.cpp
  src/SearchCubes.cpp
  src/Datastructures.cpp
  src/ObjectRegistry.cpp
  src/cgal/CGALHelper.cpp
  )


target_link_libraries(Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(Core PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
