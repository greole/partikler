
target_sources(Integrator
  PRIVATE
    STLPosIntegrator.cpp
  PUBLIC
    STLPosIntegrator.hpp
    )

target_include_directories(Integrator
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
