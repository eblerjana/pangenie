set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles
	${PROGRAM_SOURCE_DIR}/emissionprobabilitycomputer.cpp
	${PROGRAM_SOURCE_DIR}/copynumber.cpp
	${PROGRAM_SOURCE_DIR}/commandlineparser.cpp
	${PROGRAM_SOURCE_DIR}/columnindexer.cpp
	${PROGRAM_SOURCE_DIR}/dnasequence.cpp
	${PROGRAM_SOURCE_DIR}/fastareader.cpp
	${PROGRAM_SOURCE_DIR}/genotypingresult.cpp
	${PROGRAM_SOURCE_DIR}/histogram.cpp
	${PROGRAM_SOURCE_DIR}/hmm.cpp
	${PROGRAM_SOURCE_DIR}/jellyfishcounter.cpp
	${PROGRAM_SOURCE_DIR}/jellyfishreader.cpp
	${PROGRAM_SOURCE_DIR}/kmerpath.cpp
	${PROGRAM_SOURCE_DIR}/pathsampler.cpp
	${PROGRAM_SOURCE_DIR}/probabilitycomputer.cpp
	${PROGRAM_SOURCE_DIR}/probabilitytable.cpp
	${PROGRAM_SOURCE_DIR}/sequenceutils.cpp
	${PROGRAM_SOURCE_DIR}/timer.cpp
	${PROGRAM_SOURCE_DIR}/transitionprobabilitycomputer.cpp
	${PROGRAM_SOURCE_DIR}/threadpool.cpp
	${PROGRAM_SOURCE_DIR}/uniquekmercomputer.cpp
	${PROGRAM_SOURCE_DIR}/uniquekmers.cpp
	${PROGRAM_SOURCE_DIR}/variant.cpp
	${PROGRAM_SOURCE_DIR}/variantreader.cpp)

add_executable(PanGenie pggtyper.cpp ${ProjectFiles})

target_link_libraries(PanGenie ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie ${JELLYFISH_LIBRARIES})
