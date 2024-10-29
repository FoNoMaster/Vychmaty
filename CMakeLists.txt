cmake_minimum_required(VERSION 3.27)
project(Vychmaty)

add_subdirectory(src)
add_subdirectory(third_party)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()