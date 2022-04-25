cmake_minimum_required(VERSION 3.1)
project(PanGenie)

set(CMAKE_CXX_STANDARD 20)
include(FetchContent)
include(GNUInstallDirs)
enable_testing()

#find_package(PkgConfig REQUIRED)
#pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)


#find_package(jellyfish 2.3.0  QUIET)
#if (NOT jellyfish_FOUND)
#  FetchContent_Declare(
#    jellyfish
#    GIT_REPOSITORY https://github.com/gmarcais/Jellyfish
#    GIT_TAG v2.3.0)
#
#  FetchContent_GetProperties(jellyfish)
#  if (NOT jellyfish_POPULATED)
#    FetchContent_Populate(jellyfish)
#    add_subdirectory(
#      ${jellyfish_SOURCE_DIR}
#      ${jellyfish_BINARY_DIR}
#      EXCLUDE_FROM_ALL)
#    add_library(jellyfish::jellyfish ALIAS jellyfish)
#  endif ()
#endif ()

set(JUST_INSTALL_CEREAL ON CACHE BOOL "skip boost")
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



link_directories(/cluster/work/pausch/alex/software/jellyfish-2.3.0/lib cereal::cereal)
#link_directories(jellyfish::jellyfish cereal::cereal)
#include_directories(${CMAKE_BINARY_DIR}/src /cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src/include)
include_directories(/cluster/work/pausch/alex/software/jellyfish-2.3.0/include ${CMAKE_BINARY_DIR}/src /cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src/include)
#add_compile_options(${JELLYFISH_CFLAGS_OTHER})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -std=gnu++2a -O3 -flto -pipe -mavx2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -std=gnu++2a -O1 -g3")
link_libraries(pthread jellyfish-2.0)
add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
