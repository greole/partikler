project(Kernel LANGUAGES CXX)

#find_package(GTest REQUIRED)
find_package(yaml-cpp REQUIRED)

add_library(Kernel SHARED
  src/stl/Wendland2D.cpp
  src/Wendland.cpp
  )

add_dependencies(Kernel Core)

target_link_libraries(Kernel Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(Kernel PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
