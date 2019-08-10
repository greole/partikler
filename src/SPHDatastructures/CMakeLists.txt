project(SPHDatastructures LANGUAGES CXX)

#include( ${CGAL_USE_FILE} )  # needed for the linker
#find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(GTest REQUIRED)

add_library(SPHDatastructures SHARED
  src/CGALHelper.cpp
  src/SearchCubes.cpp
  src/SPHDatastructures.cpp
  src/SPHObject.cpp
  src/FileIO.cpp
  src/Helper.cpp
  )

target_include_directories(SPHDatastructures PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests)
