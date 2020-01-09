project(IO LANGUAGES CXX)


IF(WITH_H5HUT)
    IF(BUILD_HDF5)
          BUILD(HDF5)

          set(HDF5_EXTERNALLY_CONFIGURED 1)

          # set(HDF5_EXPORTED_TARGETS "HDF5")
          set(HDF5_LIBRARIES_TO_EXPORT "hdf5")

          add_subdirectory(
            ${CMAKE_CURRENT_BINARY_DIR}/HDF5-src # source dir
            ${CMAKE_CURRENT_BINARY_DIR}/HDF5-build
            )



          # EXTERNALPROJECT_GET_PROPERTY(HDF5-src install_dir)
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

    BUILD(H5Hut)

    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/H5hut-src/configure
        COMMAND ./autogen.sh
        COMMAND ./configure
                --with-hdf5=${HDF5_LIBRARY_DIRS}/..
                --prefix=${CMAKE_CURRENT_BINARY_DIR}
                --enable-parallel
        COMMAND make
        COMMAND make install
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/H5Hut-src"
        DEPENDS ${H5HutDep}
    )

    ADD_CUSTOM_TARGET(H5HUT
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/H5hut-src/configure 
      )

    add_dependencies(H5HUT hdf5_cpp-shared)

    #     EXTERNALPROJECT_GET_PROPERTY(H5HUT-download install_dir)
    #     INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS} ${H5Hut_INCLUDE_DIR})
    #     ADD_LIBRARY(H5Hut STATIC IMPORTED)
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


add_library(IO
  src/ParticleGenerator.cpp # dummy for CMAKE V3.10
  )

#TODO In CMAKE > 3.10 header path can be simplified
target_sources(IO
  PRIVATE
  "src/ParticleGenerator.cpp"
  "src/WriterBase.cpp"
  "src/SuperSPHWriter.cpp"
  "src/HDF5Writer.cpp"
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleGenerator.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/WriterBase.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/SuperSPHWriter.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HDF5Writer.hpp"
  )

MESSAGE("CMAKE CURRENT SOURCE ${CMAKE_CURRENT_SOURCE_DIR}")

add_dependencies(IO Core H5HUT)

target_link_libraries(IO
  PUBLIC
    Core
    ${YAML_CPP_LIBRARIES}
  PRIVATE
    OpenMP::OpenMP_CXX
    hdf5
    H5hut
  )

target_link_directories(IO
  PUBLIC
  "${CMAKE_CURRENT_BINARY_DIR}/lib"
  )

target_include_directories(IO
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )


#add_subdirectory(tests/SPHCore)
