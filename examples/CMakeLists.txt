cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET example)

file(GLOB TARGET_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)
add_executable(${TARGET} ${TARGET_SRC})

target_include_directories(${TARGET}
    PRIVATE
        ${GSL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/qengine/include
        ${CMAKE_SOURCE_DIR}/qengine/utils
)

target_compile_features(${TARGET} PUBLIC cxx_std_17)