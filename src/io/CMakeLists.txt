project(IO LANGUAGES CXX)

# find_package(GTest REQUIRED)
# find_package(yaml-cpp REQUIRED)

add_library(IO SHARED
  src/ParticleGenerator.cpp
  src/WriterBase.cpp
  src/SuperSPHWriter.cpp
  )

add_dependencies(IO Core)


target_link_libraries(IO Core OpenMP::OpenMP_CXX ${YAML_CPP_LIBRARIES})
target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR}/src)

#add_subdirectory(tests/SPHCore)
