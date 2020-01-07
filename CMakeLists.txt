cmake_minimum_required(VERSION 3.10)


########################
# Compiler settings   #
########################

project(PARTIKLER LANGUAGES CXX)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

MESSAGE("BUILD TYPE ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb -fPIC ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
  -Wall -Wpedantic -Wextra -march=native -fopenmp -fPIC"
  )

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
      -fno-omit-frame-pointer \
      -fsanitize=address"
    )

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
    -fno-omit-frame-pointer \
    -Rpass=loop-vectorize \
      -ffast-math \
      "
      #-Rpass-analysis=loop-vectorize \
      #-Rpass-missed=loop-vectorize \
      )

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

    add_definitions(-DWITH_GNU)

    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
      --coverage\
      -fprofile-arcs \
      -fsanitize=address \
      -ftest-coverage"
)

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
        -ffast-math \
        -march=native  \
        -ftree-vectorize \
        -g -rdynamic \
        "
        # -fopt-info-vec-missed"
        )
endif()


########################
# 3rd Party Libraries  #
########################

set(CMAKE_ARGS "${CMAKE_ARGS} \
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} \
        -DYAP_BUILD_TESTS=OFF \
        -DYAP_BUILD_PERF=OFF  \
        -DYAP_BUILD_DOC=OFF   \
        -DYAML_CPP_BUILD_TESTS=OFF \
        -DYAML_CPP_BUILD_TOOLS=OFF \
        -DYAML_CPP_BUILD_CONTRIB=OFF "
  )

find_package(OpenMP REQUIRED)

include(ExternalProject)

include(cmake/build_3rd_party.cmake)

OPTION(BUILD_CGAL "Building CGAL" OFF)
IF(BUILD_CGAL)
    BUILD(CGAL)
ENDIF()

find_package(CGAL CONFIG REQUIRED)

OPTION(BUILD_GOOGLE_TEST "Building google test" OFF)

# MESSAGE("ADDING GOOGLE TEST")
# add_subdirectory(3rdParty/googletest)
# MESSAGE("ADDING GOOGLE BENCHMARK")
enable_testing()
IF(NOT BUILD_GOOGLE_TEST)
    MESSAGE("Searching for GTest")
    FIND_PACKAGE(GTest REQUIRED)
ENDIF()


IF(BUILD_GOOGLE_TEST)
    BUILD(googletest)

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # # Add googletest directly to our build. This defines
    # # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)

    #  message("BUILDING GTEST DONE")

    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_BINARY_DIR}/googletest-build/lib/cmake/GTest/")
ENDIF()

find_package(GTest)


OPTION(BUILD_YAML_CPP "Building yaml-cpp parser" OFF)

IF(NOT BUILD_YAML_CPP)
    MESSAGE("Searching for YAML-CPP")
    FIND_PACKAGE(YAML-CPP CONFIG REQUIRED)
    # SET(YAML_CPP_LIBRARY yaml-cpp)
    MESSAGE("Found YAML-CPP CONFIG UNDER" ${YAML-CPP_DIR})
    MESSAGE("Found " ${YAML_CPP_LIBRARIES})
ENDIF()

IF(BUILD_YAML_CPP)

    SET(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    SET(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    SET(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
    BUILD(YAML)
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-src
                     ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-build
                     EXCLUDE_FROM_ALL)
ENDIF()


OPTION(BUILD_YAP "Building YAP" ON)
IF(NOT BUILD_YAP)
    MESSAGE("Searching for YAP")
    FIND_PACKAGE(yap REQUIRED)
ENDIF()

IF(BUILD_YAP)
    SET(YAP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    SET(YAP_BUILD_PERF OFF CACHE BOOL "" FORCE)
    SET(YAP_BUILD_DOC OFF CACHE BOOL  "" FORCE)
    BUILD(YAP)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    # add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/yap-src
    #                  ${CMAKE_CURRENT_BINARY_DIR}/yap-build
    #                  EXCLUDE_FROM_ALL)

    include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/yap-src/include)
    message("BUILDING YAP DONE")
ENDIF()


OPTION(WITH_H5HUT "With H5Hut" ON)
OPTION(BUILD_HDF5 "Building hdf5" OFF)
OPTION(BUILD_H5HUT "Build H5Part data-Format" ON)
IF(WITH_H5HUT)
    IF(BUILD_HDF5)
          BUILD(HDF5)

          add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/HDF5-src
                          ${CMAKE_CURRENT_BINARY_DIR}/HDF5-build
                          EXCLUDE_FROM_ALL)
          # EXTERNALPROJECT_GET_PROPERTY(HDF5-download install_dir)
          # INCLUDE_DIRECTORIES(SYSTEM ${install_dir}/include)
          # ADD_LIBRARY(HDF5 STATIC IMPORTED)
          # SET_TARGET_PROPERTIES(HDF5
          #     PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libhdf5.a)
          # SET(HDF5_LIBRARIES hdf5 dl z m)
          # SET(HDF5_LIBRARY_DIRS  ${CMAKE_BINARY_DIR})
      ELSE()
          LIST(APPEND CMAKE_PREFIX_PATH "$ENV{LD_LIBRARY_PATH}")
          FIND_PACKAGE(HDF5)
          # IF(NOT HDF5_FOUND)
          #     MESSAGE("No HDF5 compiled with cmake found")
          #     FIND_PATH(HDF5_LIBRARY_DIRS libhdf5.a
          #         HINTS "/usr/lib/x86_64-linux-gnu/hdf5/openmpi" "/usr/lib/" $ENV{LD_LIBRARY_PATH}
          #       )
          #   MESSAGE(${HDF5_LIBRARY_DIRS})
          #   INCLUDE_DIRECTORIES("${HDF5_LIBRARY_DIRS}/../include")
          #   LINK_DIRECTORIES("${HDF5_LIBRARY_DIRS}/../lib")
          #   SET(HDF5_LIBRARIES hdf5 dl z m)
          # ENDIF()
      SET(H5HutDep "")
    ENDIF()

    IF(NOT BUILD_H5HUT)
        FIND_PACKAGE(H5HUT REQUIRED)
    ENDIF()

    # IF(NOT H5HUT_FOUND)
    #     MESSAGE("Build H5Hut")
    #     BUILD(H5Hut)

    #     ADD_CUSTOM_COMMAND(
    #         OUTPUT ${PROJECT_SOURCE_DIR}/ThirdParty/H5hut/configure
    #         COMMAND ./autogen.sh
    #         COMMAND ./configure --with-hdf5=${HDF5_LIBRARY_DIRS}/.. --prefix=${CMAKE_BINARY_DIR} --enable-parallel && make && make install
    #         WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/ThirdParty/H5hut/"
    #         DEPENDS ${H5HutDep}
    #     )

    #     ADD_CUSTOM_TARGET(H5HUTBUILD
    #         DEPENDS ${PROJECT_SOURCE_DIR}/ThirdParty/H5hut/configure)

    #     EXTERNALPROJECT_GET_PROPERTY(H5HUT-download install_dir)
    #     INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS} ${H5Hut_INCLUDE_DIR})
    #     ADD_LIBRARY(H5Hut STATIC IMPORTED)
    #     ADD_DEPENDENCIES(H5Hut H5HUTBUILD)
    #     SET_TARGET_PROPERTIES(H5Hut
    #         PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libH5hut.a)
    #     SET(H5Hut_LIBRARY H5Hut)
    # ENDIF()

    # SET(HDF5_PREFER_PARALLEL TRUE)
    # ADD_COMPILE_DEFINITIONS(WITH_H5HUT)
    # IF (HAVE_API2_FUNCTIONS)
    #     add_compile_definitions(WITH_H5HUT_V2)
    # ENDIF (HAVE_API2_FUNCTIONS)
ENDIF()


########################
# gcov                 #
########################

IF ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
include (cmake/CodeCoverage.cmake)

append_coverage_compiler_flags()

setup_target_for_coverage_lcov(
  NAME coverage
  EXECUTABLE ./src/core/tests/Field/Partikler_Test_Field
  BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src"
  EXCLUDE "build/*" "/usr/*"
  )

# setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
# target_link_libraries(${PROJECT_TEST_NAME} gcov)

ENDIF()
########################
# ctags                #
########################


OPTION(GENERATE_TAGS "Generating tags" OFF)
IF(GENERATE_TAGS)
    add_custom_command ( OUTPUT tags
        COMMAND ctags -R ${CMAKE_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
ELSE()
    add_custom_command ( OUTPUT tags
        COMMAND touch tags
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDIF()

########################
# own Libraries        #
########################

add_subdirectory(src/core)
add_subdirectory(src/boundary_conditions)
add_subdirectory(src/neighbours)
add_subdirectory(src/io)
add_subdirectory(src/kernel)
add_subdirectory(src/transporteqn)
add_subdirectory(src/integrator)
add_subdirectory(src/geometry)

add_executable(partikler src/partikler.cpp tags)

#target_link_libraries(partikler SPHDatastructures ${YAML_CPP_LIBRARIES})

# target_link_libraries(partikler BoundaryConditions Geometry Integrator Neighbours Core TransportEqn Kernel IO OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_link_libraries(partikler
  Core
  Neighbours
  BoundaryConditions
  IO
  TransportEqn
  Kernel
  Integrator
  Geometry
  OpenMP::OpenMP_CXX
  ${YAML_CPP_LIBRARIES}
  ${HDF5_LIBRARIES}
  )


OPTION(COVERALLS "Calling Coveralls" OFF)

IF(COVERALLS AND ${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" AND ${CMAKE_BUILD_TYPE} MATCHES "Debug")
      add_custom_command ( OUTPUT covered
          COMMAND coveralls && touch covered
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
ELSE()
    add_custom_command ( OUTPUT covered
        COMMAND touch covered
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDIF()

# add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
# add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)
#
# # Integeration tests
# IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#   add_test(valgrind valgrind --tool=callgrind --include=.. ./partikler ../example/bunny.stl bunny 0.5)
# ELSE()
#   add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)
# ENDIF()
