cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET unit_tests)
project(${TARGET})

include(cmake/functions.cmake)

add_subdirectory(qengine/utils)
add_subdirectory(qengine/include)