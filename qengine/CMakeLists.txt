cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET qengine)
project(${TARGET})

set(QENGINE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Folder for Quantum Engine")

include_directories(${TARGET}
    PRIVATE
        ${GSL_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/utils
)

file(GLOB TARGET_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(${TARGET} ${TARGET_SRC})

target_compile_features(${TARGET} PUBLIC cxx_std_17)
