add_library(PGGTyperLib SHARED emissionprobabilitycomputer.cpp copynumber.cpp)

add_executable(PGGTyper main.cpp)

target_link_libraries(PGGTyper PGGTyperLib)
