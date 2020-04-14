
target_sources(Core
  PRIVATE
    CGALSearchCubes.cpp
    CGALHelper.cpp
    CGALField.cpp
  PUBLIC
    CGALSearchCubes.hpp
    CGALHelper.hpp
    CGALTYPEDEFS.hpp
    CGALField.hpp
    )

target_link_libraries(Core
  PUBLIC
    ${CGAL_LIBRARIES}
  )

target_include_directories(Core PUBLIC ${CMAKE_CURRENT_LIST_DIR})
