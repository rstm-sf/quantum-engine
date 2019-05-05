cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET example)

file(GLOB TARGET_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)
add_executable(${TARGET} ${TARGET_SRC})

target_compile_features(${TARGET} PUBLIC cxx_std_17)

target_link_libraries(${TARGET}
    PRIVATE
        GSL
        qengine
)
