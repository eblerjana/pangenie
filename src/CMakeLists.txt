add_library(PanGenieLib SHARED 
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	commandlineparser.cpp
	columnindexer.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	histogram.cpp
	hmm.cpp
	jellyfishcounter.cpp
	jellyfishreader.cpp
	kmerpath.cpp
	kmerparser.cpp
	pathsampler.cpp
	probabilitycomputer.cpp
	probabilitytable.cpp
	sequenceutils.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	uniquekmers.cpp
	variant.cpp
	variantreader.cpp)

add_executable(PanGenie pangenie.cpp)

target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LIBRARIES})