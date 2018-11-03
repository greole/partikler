cmake_minimum_required(VERSION 3.10)
project( PARTIKLER LANGUAGES CXX)

find_package(CGAL  REQUIRED)
include( ${CGAL_USE_FILE} )  # needed for the linker

find_package(GTest REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE "RELEASE")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wpedantic -O3 -march=native")

enable_testing()
include_directories(include)

add_subdirectory(src/SPHDatastructures)
add_subdirectory(src/SPHio)
add_subdirectory(src/SPHCore)

add_executable(partikler src/partikler.cpp)

target_link_libraries(partikler SPHDatastructures SPHIO SPHCore)

add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)

# Integeration tests
IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  add_test(valgrind valgrind --tool=callgrind --include=.. ./partikler ../example/bunny.stl bunny 0.5)
ELSE()
  add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)
ENDIF()



# include(CTest)
#
# set(PROJECT_SRC_DIR ${PROJECT_SOURCE_DIR}/src)
#
# include_directories (BEFORE "${PROJECT_SRC_DIR}")
#
# add_subdirectory(tests)
#
#
#
# #set(CMAKE_CXX_COMPILER "clang++")
# #set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -O3  -march=native") #-stdlib=libc++ -std=c++1z")
#
#
#
# OPTION (USE_OpenMP "Use OpenMP" ON)
# IF(USE_OpenMP)
#   FIND_PACKAGE(OpenMP)
#   IF(OPENMP_FOUND)
#     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#   ENDIF()
# ENDIF()
#
# # set(CMAKE_CXX_FLAGS "-Wall -O3  -march=native ") #-stdlib=libc++ -std=c++1z")
#
# if ( CGAL_FOUND )
#
#   include( ${CGAL_USE_FILE} )
#
#   # find_package(TBB QUIET)
#   # include(${TBB_USE_FILE})
#   list(APPEND CGAL_3RD_PARTY_LIBRARIES)
#
#   include( CGAL_CreateSingleSourceCGALProgram )
#
#   include(src/git_watcher.cmake)
#
#   FILE(GLOB src src/io/*.h src/*.h src/*.cpp )
#
#   add_executable(partikler partikler.cpp ${SRC})
#
#   add_dependencies(partikler AlwaysCheckGit)
#
#   # Unit tests
# else()
#
#     message(STATUS "This program requires the CGAL library, and will not be compiled.")
#
# endif()
#
