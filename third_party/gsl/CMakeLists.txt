configure_file(CMakeLists.txt.in gsl-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gsl-download )
if(result)
  message(FATAL_ERROR "CMake step for GSL failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gsl-download )
if(result)
  message(FATAL_ERROR "Build step for GSL failed: ${result}")
endif()

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/gsl-src
                 ${CMAKE_CURRENT_BINARY_DIR}/gsl-build)

set(GSL_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsl-src/include CACHE INTERNAL "Include folder for GSL")