FUNCTION(BUILD NAME)
    configure_file(cmake/${NAME}.in ${NAME}-download/CMakeLists.txt)

    message("CMAKE ${NAME}")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${NAME}-download )

    if(result)
      message(FATAL_ERROR "CMake step for ${NAME} failed: ${result}")
    endif()

    message("BUILDING ${NAME}")
    execute_process(
      COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${NAME}-download )
    if(result)
      message(FATAL_ERROR "Build step for ${NAME} failed: ${result}")
    endif()

ENDFUNCTION(BUILD)
