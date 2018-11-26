cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET include_tests)

custom_add_executable_from_dir(${TARGET})

target_include_directories(${TARGET}
    PRIVATE
        ${GSL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/qengine/include
        ${CMAKE_SOURCE_DIR}/qengine/utils
)

target_compile_features(${TARGET} PUBLIC cxx_std_17)

add_test(test_all ${TARGET})