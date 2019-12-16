cmake_minimum_required(VERSION 3.10)


########################
# Compiler settings   #
########################

project(PARTIKLER LANGUAGES CXX)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "RELEASE")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
    -Wall -Wpedantic -Wextra -march=native -fopenmp -fPIC"
    )

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
      -fno-omit-frame-pointer \
      -Rpass=loop-vectorize \
      -ffast-math \
      -fsanitize=address \
      "
      #-Rpass-analysis=loop-vectorize \
      #-Rpass-missed=loop-vectorize \
      )

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

    add_compile_definitions(WITH_GNU)

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
        -ffast-math \
        -march=native  \
        -ftree-vectorize \
        "
        # -fopt-info-vec-missed"
        )
endif()


########################
# 3rd Party Libraries  #
########################
#add_subdirectory(3rdParty/cgal)

include(ExternalProject)

OPTION(BUILD_CGAL "Building CGAL" OFF)
IF(BUILD_CGAL)
    message("BUILDING CGAL")
    # Download and unpack googletest at configure time
    configure_file(cmake/CGAL.in cgal-download/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cgal-download )
    if(result)
      message(FATAL_ERROR "CMake step for cgal failed: ${result}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cgal-download )
    if(result)
      message(FATAL_ERROR "Build step for cgal failed: ${result}")
    endif()

    execute_process(COMMAND mkdir -p build
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cgal-src )
    if(result)
      message(FATAL_ERROR "Build step for cgal failed: ${result}")
    endif()

    execute_process(COMMAND cmake  -DCGAL_HEADER_ONLY=ON -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/cgal-src/build ..
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cgal-src/build )
    if(result)
      message(FATAL_ERROR "Build step for cgal failed: ${result}")
    endif()

    execute_process(COMMAND make install
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cgal-src/build )
    if(result)
      message(FATAL_ERROR "Build step for cgal failed: ${result}")
    endif()

    # add_compile_definitions(CGAL_HEADER_ONLY)
    #SET_TARGET_PROPERTIES(CGAL
    #    PROPERTIES
    #    CGAL_HEADER_ONLY TRUE)

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR}/cgal-src)
    set(CGAL_DIR  ${CMAKE_CURRENT_BINARY_DIR}/cgal-src)
    find_package(CGAL  REQUIRED)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/cgal-src/build/include)
ENDIF()


# find_package(CGAL  REQUIRED)

# MESSAGE("Found CGAL CONFIG UNDER" ${CGAL_DIR})
# include( ${CGAL_USE_FILE} )  # needed for the linker


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
    # Download and unpack googletest at configure time
    configure_file(cmake/GTest.in googletest-download/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
      message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()

    message("BUILDING GTEST")
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
      message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)

     message("BUILDING GTEST DONE")

    execute_process(COMMAND make
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-src
      )
    find_package(GTest)
ENDIF()


#FIND_PACKAGE(GTest REQUIRED)

# add_subdirectory(3rdParty/benchmark)

# build tests (targets: gtest_main, gtest)
#set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)

INCLUDE(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
OPTION(BUILD_YAML_CPP "Building yaml-cpp parser" OFF)

IF(NOT BUILD_YAML_CPP)
    MESSAGE("Searching for YAML-CPP")
    FIND_PACKAGE(YAML-CPP CONFIG REQUIRED)
    # SET(YAML_CPP_LIBRARY yaml-cpp)
    MESSAGE("Found YAML-CPP CONFIG UNDER" ${YAML-CPP_DIR})
    MESSAGE("Found " ${YAML_CPP_LIBRARIES})
ENDIF()

IF(BUILD_YAML_CPP)
    MESSAGE("BUILDING YAML_CPP ")
      EXTERNALPROJECT_ADD(yaml-cpp-download
          GIT_REPOSITORY    "https://github.com/jbeder/yaml-cpp"
          # GIT_TAG           "0.6.2"
          DOWNLOAD_DIR      "${PROJECT_SOURCE_DIR}/3rdParty/"
          PREFIX            "${CMAKE_BINARY_DIR}"
          SOURCE_DIR        "${PROJECT_SOURCE_DIR}/3rdParty/yaml-cpp"
          CMAKE_ARGS        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}
                            -DYAML_CPP_BUILD_TESTS=OFF
                            -DYAML_CPP_BUILD_TOOLS=OFF
                            -DYAML_CPP_BUILD_CONTRIB=OFF
                            -DYAML_CPP_INSTALL=ON
                            -DBUILD_SHARED_LIBS=ON
        )

      EXTERNALPROJECT_GET_PROPERTY(yaml-cpp-download install_dir)
      INCLUDE_DIRECTORIES(${install_dir}/include)

      ADD_LIBRARY(yaml-cpp SHARED IMPORTED)
      SET_TARGET_PROPERTIES(yaml-cpp
          PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libyaml-cpp.so)
      SET(YAML_CPP_LIBRARY yaml-cpp)
      MESSAGE("SET YAML_CPP_LIBRARY " ${YAML_CPP_LIBRARY})
ENDIF()


OPTION(BUILD_YAP "Building YAP" ON)
# IF(NOT BUILD_YAP)
#     MESSAGE("Searching for YAP")
#     FIND_PACKAGE(yap REQUIRED)
# ENDIF()

IF(BUILD_YAP)
    # Download and unpack googletest at configure time
    configure_file(cmake/YAP.in yap-download/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/yap-download )
    if(result)
      message(FATAL_ERROR "CMake step for yap failed: ${result}")
    endif()

    message("BUILDING YAP")
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/yap-download )
    if(result)
      message(FATAL_ERROR "Build step for yap failed: ${result}")
    endif()

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/yap-src
                     ${CMAKE_CURRENT_BINARY_DIR}/yap-build
                     EXCLUDE_FROM_ALL)

                 message("BUILDING YAP DONE")

    execute_process(COMMAND make
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/yap-src
      )

    execute_process(COMMAND make
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/yap-src
      )

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/yap-src/include)
ENDIF()

########################
# own Libraries        #
########################

add_subdirectory(src/core)
add_subdirectory(src/boundary_conditions)
add_subdirectory(src/neighbours)
# add_subdirectory(src/SPHDatastructures)
# add_subdirectory(src/kernel)
# add_subdirectory(src/transporteqn)
# add_subdirectory(src/integrator)
# add_subdirectory(src/geometry)
# add_subdirectory(src/io)

add_executable(partikler src/partikler.cpp tags)

#target_link_libraries(partikler SPHDatastructures ${YAML_CPP_LIBRARIES})

# target_link_libraries(partikler BoundaryConditions Geometry Integrator Neighbours Core TransportEqn Kernel IO OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_link_libraries(partikler Core Neighbours BoundaryConditions OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})



# add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
# add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)
#
# # Integeration tests
# IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#   add_test(valgrind valgrind --tool=callgrind --include=.. ./partikler ../example/bunny.stl bunny 0.5)
# ELSE()
#   add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)
# ENDIF()
