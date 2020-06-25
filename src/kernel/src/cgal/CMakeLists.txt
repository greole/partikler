
target_sources(Kernel
  PRIVATE
    Wendland2D.cpp
  PUBLIC
    Wendland2D.hpp
    )

target_include_directories(Kernel
  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
