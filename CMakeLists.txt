cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(quantum-engine)

add_subdirectory(qengine)

#############################################
# Unit tests

add_subdirectory(third_party/googletest)

# enable CTest testing
enable_testing()
set(BUILD_TESTING ON)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()