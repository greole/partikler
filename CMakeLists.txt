# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
project( Partikler )

include(CTest)

cmake_minimum_required(VERSION 2.8.10)

include(git_watcher.cmake)

find_package(CGAL QUIET)


set(CMAKE_CXX_FLAGS "-Wall -O3 -fopenmp  -march=native -static-libstdc++ -D_GLIBCXX_PARALLEL") #-stdlib=libc++ -std=c++1z")
# set(CMAKE_CXX_FLAGS "-Wall -O3  -march=native ") #-stdlib=libc++ -std=c++1z")

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  # find_package(TBB QUIET)
  # include(${TBB_USE_FILE})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES)

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  add_executable(partikler particle_helper.h io.h git.h partikler.cpp)


  add_dependencies(partikler AlwaysCheckGit)

  add_test(BuildBunny1   time ./partikler ../example/bunny.stl bunny 1.0)
  add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
  add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

