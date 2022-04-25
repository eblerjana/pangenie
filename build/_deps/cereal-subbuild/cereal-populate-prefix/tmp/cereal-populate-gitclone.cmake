
if(NOT "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt" IS_NEWER_THAN "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/cluster/apps/sfos/bin/git"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/USCiLab/cereal" "cereal-src"
    WORKING_DIRECTORY "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/USCiLab/cereal'")
endif()

execute_process(
  COMMAND "/cluster/apps/sfos/bin/git"  checkout v1.3.1 --
  WORKING_DIRECTORY "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'v1.3.1'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/cluster/apps/sfos/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt"
    "/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/cluster/work/pausch/alex/software/Alex-pangenie/build/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitclone-lastrun.txt'")
endif()

