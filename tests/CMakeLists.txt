add_executable(tests tests.cpp EmissionProbabilityComputerTest.cpp CopyNumberTest.cpp)

target_link_libraries(tests PGGTyperLib)
