cmake_minimum_required(VERSION 3.1)
project(PanGenie)

set(CMAKE_CXX_STANDARD 20)
include(FetchContent)
include(GNUInstallDirs)
enable_testing()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)

find_package(cereal 1.3.1 QUIET)
if (NOT cereal_FOUND)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal
    GIT_TAG v1.3.1)

  FetchContent_GetProperties(cereal)
  if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(
      ${cereal_SOURCE_DIR}
      ${cereal_BINARY_DIR}
      EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()



link_directories(${JELLYFISH_LIBRARY_DIRS} cereal:cereal)

include_directories(${JELLYFISH_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/src /cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src/include)
#add_compile_options(${JELLYFISH_CFLAGS_OTHER})
#target_link_libraries(PanGenie
#  cereal::cereal)

#target_include_directories(PanGenie PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -std=gnu++2a -O3 -flto -pipe -mavx2")

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
