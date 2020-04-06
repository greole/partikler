project(Core LANGUAGES CXX)

find_package(yaml-cpp REQUIRED)

# TODO move SearchCubes to neighbours?

add_library(Core OBJECT)

target_sources(Core
  PRIVATE
    src/Scalar.cpp
    src/Vec3.cpp
    src/Object.cpp
    src/Field.cpp
    src/FieldOps.cpp
    src/Models.cpp
    src/Helper.cpp
    src/SearchCubes.cpp
    src/Time.cpp
    src/Equation.cpp
    src/ObjectRegistry.cpp
  PUBLIC
    src/Scalar.hpp
    src/Vec3.hpp
    src/Object.hpp
    src/Field.hpp
    src/FieldOps.hpp
    src/Models.hpp
    src/Helper.hpp
    src/Time.hpp
    src/Equation.hpp
    src/SearchCubes.hpp
    src/ObjectRegistry.hpp
  )


target_compile_definitions(Core
  PUBLIC
  SINGLEPRECISION=${SINGLEPRECISIONVAL}
  SINGLEPRECISION=${DOUBLEPRECISIONVAL}
  )

target_link_libraries(Core
  PUBLIC
    ${YAML_CPP_LIBRARIES}
    OpenMP::OpenMP_CXX
  )

target_include_directories(Core
  PUBLIC
    src
  )

IF(CGAL_SUPPORT)
  add_subdirectory(src/cgal)
ENDIF()

add_subdirectory(tests)

