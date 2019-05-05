cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET qengine)
project(${TARGET} CXX)

include(ExternalProject)
find_package(Git)

#set(QENGINE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Folder for Quantum Engine")

# creates a library Quantum Engine which is an interface (header files only)
add_library(${TARGET} INTERFACE)

# determine whether this is a standalone project or included by other projects
set(QENGINE_STANDALONE_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(QENGINE_STANDALONE_PROJECT ON)
endif ()

set(QENGINE_CXX_STANDARD "14" CACHE STRING "Use c++ standard")
set(QENGINE_CXX_STD "cxx_std_${QENGINE_CXX_STANDARD}")

# when minimum version required is 3.8.0 remove if below
# both branches do exactly the same thing
target_compile_features(GSL INTERFACE "${QENGINE_CXX_STD}")
# on *nix systems force the use of -std=c++XX instead of -std=gnu++XX (default)
set(CMAKE_CXX_EXTENSIONS OFF)

# add definitions to the library and targets that consume it
target_compile_definitions(${TARGET} INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:
        # remove unnecessary warnings about unchecked iterators
        _SCL_SECURE_NO_WARNINGS
        # remove deprecation warnings about std::uncaught_exception() (from catch)
        _SILENCE_CXX17_UNCAUGHT_EXCEPTION_DEPRECATION_WARNING
    >
)

# add include folders to the library and targets that consume it
target_include_directories(${TARGET} INTERFACE
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
    >
)
