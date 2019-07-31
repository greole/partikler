project(SPHDatastructures LANGUAGES CXX)

find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

add_library(SPHDatastructures SHARED
  src/SPHDatastructures.cpp
  src/SPHObject.cpp
  src/SearchCubes.cpp
  src/FileIO.cpp
  src/CGALHelper.cpp
  src/Helper.cpp
  )

target_include_directories(SPHDatastructures PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(tests)
