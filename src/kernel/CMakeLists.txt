project(Kernel LANGUAGES CXX)

#find_package(GTest REQUIRED)
find_package(yaml-cpp REQUIRED)

add_library(Kernel SHARED)

target_sources(Kernel
  PRIVATE
  src/Wendland.cpp
  src/Gauss.cpp
  PUBLIC
  src/Wendland.hpp
  src/Gauss.hpp
  )


add_dependencies(Kernel Core)

target_link_libraries(Kernel
  PUBLIC
  Core
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES}
)

target_include_directories(Kernel
  PUBLIC
  src
  )

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

#add_subdirectory(tests/SPHCore)
