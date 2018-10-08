# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Generator_Examples )

cmake_minimum_required(VERSION 2.8.10)

include(git_watcher.cmake)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  find_package(TBB QUIET)

  include(${TBB_USE_FILE})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  add_executable(partikler  git.h partikler.cpp)

  add_dependencies(partikler AlwaysCheckGit)

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

