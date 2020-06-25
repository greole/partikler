
target_sources(BoundaryConditions
  PRIVATE
    GenerateBoundaryParticles.cpp
  PUBLIC
   GenerateBoundaryParticles.hpp
  )

target_include_directories(BoundaryConditions
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
