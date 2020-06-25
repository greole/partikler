project(IO LANGUAGES CXX)

IF(WITH_H5HUT)
    IF(BUILD_HDF5)

          EXTERNALPROJECT_ADD(HDF5-download
            DOWNLOAD_NAME     "hdf5-1-10-5.tar.gz"
            URL_HASH          "MD5=7c19d6b81ee2a3ba7d36f6922b2f90d3"
            URL               "https://www.hdfgroup.org/package/hdf5-1-10-5-tar-bz2/?wpdmdl=13570&refresh=5c88e0ec5645e1552474348%"
            # PREFIX            "${CMAKE_CURRENT_BINARY_DIR}"
            SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/HDF5-src"
            BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/HDF5-build"
            CMAKE_ARGS        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
                              -DBUILD_TESTING=FALSE
                              -DHDF5_BUILD_TOOLS=FALSE
                              -DHDF5_BUILD_EXAMPLES=FALSE
                              #-DALLOW_UNSUPPORTED=TRUE
                              #-DHDF5_ENABLE_PARALLEL=TRUE
          )

      ELSE()
          LIST(APPEND CMAKE_PREFIX_PATH "$ENV{LD_LIBRARY_PATH}")
          FIND_PACKAGE(HDF5)
    ENDIF()

    IF(NOT BUILD_H5HUT)
        FIND_PACKAGE(H5HUT REQUIRED)
    ENDIF()

    BUILD(H5Hut)

    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/H5hut-src/configured
        COMMAND ./autogen.sh
        COMMAND ./configure
                CPPFLAGS=-fPIC
                --with-hdf5=${CMAKE_CURRENT_BINARY_DIR}
                --prefix=${CMAKE_CURRENT_BINARY_DIR}
                #--enable-shared
                #--enable-parallel
        COMMAND make
        COMMAND make install
        COMMAND touch configured
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/H5Hut-src"
    )

    ADD_CUSTOM_TARGET(H5HUT
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/H5hut-src/configured
      )

    add_dependencies(H5HUT HDF5-download)
ENDIF()

add_library(IO SHARED)

target_sources(IO
  PRIVATE
  src/WriterBase.cpp
  src/ReaderBase.cpp
  src/SuperSPHWriter.cpp
  src/HDF5Writer.cpp
  PUBLIC
  src/WriterBase.hpp
  src/ReaderBase.hpp
  src/SuperSPHWriter.hpp
  src/HDF5Writer.hpp
  )

MESSAGE("CMAKE CURRENT SOURCE ${CMAKE_CURRENT_SOURCE_DIR}")

add_dependencies(IO Core H5HUT)

target_link_libraries(IO
  PUBLIC
    Core
    ${YAML_CPP_LIBRARIES}
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
  "${CMAKE_CURRENT_BINARY_DIR}/include"
  )

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

#add_subdirectory(tests/SPHCore)
