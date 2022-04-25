# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.20.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(cereal-populate NONE)

# Pass through things we've already detected in the main project to avoid
# paying the cost of redetecting them again in ExternalProject_Add()
set(GIT_EXECUTABLE [==[/cluster/apps/sfos/bin/git]==])
set(GIT_VERSION_STRING [==[2.15.1]==])
set_property(GLOBAL PROPERTY _CMAKE_FindGit_GIT_EXECUTABLE_VERSION
  [==[/cluster/apps/sfos/bin/git;2.15.1]==]
)


include(ExternalProject)
ExternalProject_Add(cereal-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/USCiLab/cereal" "GIT_TAG" "v1.3.1"
                    SOURCE_DIR          "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src"
                    BINARY_DIR          "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


