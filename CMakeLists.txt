cmake_minimum_required(VERSION 3.10)
project( PARTIKLER LANGUAGES CXX)

find_package(CGAL  REQUIRED)
include( ${CGAL_USE_FILE} )  # needed for the linker

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_BUILD_TYPE "RELEASE")
set(CMAKE_BUILD_TYPE "DEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wpedantic -O3 -march=native")

add_custom_command ( OUTPUT tags
    COMMAND ctags -R ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

enable_testing()

add_subdirectory(src/SPHDatastructures)
# add_subdirectory(src/SPHio)
# add_subdirectory(src/SPHCore)

add_executable(partikler src/partikler.cpp tags)

target_link_libraries(partikler SPHDatastructures) # SPHIO SPHCore)

# add_test(BuildBunny0.5 time ./partikler ../example/bunny.stl bunny 0.5)
# add_test(BuildBunny0.25 time ./partikler ../example/bunny.stl bunny 0.25)
#
# # Integeration tests
# IF (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#   add_test(valgrind valgrind --tool=callgrind --include=.. ./partikler ../example/bunny.stl bunny 0.5)
# ELSE()
#   add_test(BuildBunny0.1 time ./partikler ../example/bunny.stl bunny 0.1)
# ENDIF()
