
target_sources(Geometry
  PRIVATE
    STLLimitedDx.cpp
  PUBLIC
    STLLimitedDx.hpp
  )

target_include_directories(Geometry
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
