# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
project( PARTIKLER )

include(CTest)
find_package(GTest REQUIRED)

set(PROJECT_SRC_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories (BEFORE "${PROJECT_SRC_DIR}")

add_subdirectory(tests)

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

set(CMAKE_CXX_FLAGS "-Wall -O3 -fopenmp  -march=native -static-libstdc++ -D_GLIBCXX_PARALLEL") #-stdlib=libc++ -std=c++1z")
# set(CMAKE_CXX_FLAGS "-Wall -O3  -march=native ") #-stdlib=libc++ -std=c++1z")

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  # find_package(TBB QUIET)
  # include(${TBB_USE_FILE})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES)

  include( CGAL_CreateSingleSourceCGALProgram )

  include(src/git_watcher.cmake)

  FILE(GLOB src src/*.h)

  add_executable(partikler partikler.cpp ${SRC})

  add_dependencies(partikler AlwaysCheckGit)

  # Unit tests

  # Integeration tests
  add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
  add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)
  add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

