
target_sources(IO
  PRIVATE
    CGALParticleGenerator.cpp
  PUBLIC
    CGALParticleGenerator.hpp
)

target_link_libraries(IO
  PUBLIC
    ${CGAL_LIBRARIES}
)

target_include_directories(IO PUBLIC ${CMAKE_CURRENT_LIST_DIR})
