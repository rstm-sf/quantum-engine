cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(TARGET utils_tests)

custom_add_executable_from_dir(${TARGET})

target_include_directories(${TARGET}
    PRIVATE
        "${CMAKE_SOURCE_DIR}/qengine/utils"
)

target_compile_features(${TARGET} PUBLIC cxx_std_14)

add_test(test_all ${TARGET})