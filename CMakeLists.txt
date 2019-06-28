cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
set(CMAKE_LINKER  "/usr/bin/ldd")
set(CMAKE_C_COMPILER "/usr/local/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++")
project(PARTIKLER LANGUAGES CXX)

find_package(CGAL  REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "RELEASE")
#set(CMAKE_BUILD_TYPE "DEBUG")

set(ENV{LD_PRELOAD} /usr/local/lib64/libCGAL.so.13)
set(CMAKE_INSTALL_RPATH "/usr/local/lib64")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
    -Wall -Wpedantic -Wextra -march=native -fopenmp"
    )

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
      -fno-omit-frame-pointer \
      -Rpass=loop-vectorize \
      -Rpass-analysis=loop-vectorize \
      -Rpass-missed=loop-vectorize \
      -ffast-math \
      -fsanitize=address \
      "
      )

    set(LINK_FLAGS
      "-rpath /usr/local/lib64/ \
       -L/usr/local/lib64/ \
      ")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
        -ffast-math \
        -march=native  \
        -ftree-vectorize \
        -fopt-info-vec-missed"
        )
endif()

# add_subdirectory(3rdParty/cgal)

include( ${CGAL_USE_FILE} )  # needed for the linker

add_custom_command ( OUTPUT tags
    COMMAND ctags -R ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

enable_testing()


add_subdirectory(src/SPHDatastructures)
add_subdirectory(src/SPHRuntime)
add_subdirectory(src/SPHModels)

# build google benchmark (target: benchmark)
# do not build tests of benchmarking lib
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
add_subdirectory(3rdParty/benchmark)

# build tests (targets: gtest_main, gtest)
add_subdirectory(3rdParty/googletest)

add_executable(partikler src/partikler.cpp tags)

target_link_libraries(partikler SPHDatastructures SPHIO SPHCore OpenMP::OpenMP_CXX) # SPHIO SPHCore)

# add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
# add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)
#
# # Integeration tests
# IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#   add_test(valgrind valgrind --tool=callgrind --include=.. ./partikler ../example/bunny.stl bunny 0.5)
# ELSE()
#   add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)
# ENDIF()
