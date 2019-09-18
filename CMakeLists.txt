cmake_minimum_required(VERSION 3.1)
project(PGGTyper)

set (CMAKE_CXX_STANDARD 11)

enable_testing()

#find_package(Boost 1.67 REQUIRED system)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)

link_directories(${JELLYFISH_LIBRARY_DIRS})

#include_directories(${JELLYFISH_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/src)
include_directories(${JELLYFISH_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/src)
add_compile_options(${JELLYFISH_CFLAGS_OTHER})

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
