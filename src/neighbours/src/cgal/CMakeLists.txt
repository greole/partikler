
target_sources(Neighbours
  PRIVATE
    STLParticleNeighbours.cpp
  PUBLIC
    STLParticleNeighbours.hpp
    )

target_include_directories(Neighbours
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
