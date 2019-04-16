cmake_minimum_required(VERSION 3.6)
project(PGGTyper)

enable_testing()

#find_package(Boost 1.36 COMPONENTS program_options iostreams unit_test_framework REQUIRED)
include_directories(${CMAKE_BINARY_DIR}/src)

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)
