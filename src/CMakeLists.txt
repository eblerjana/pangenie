add_library(PGGTyperLib SHARED 
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	columnindexer.cpp
	copynumberassignment.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	hmm.cpp
	kmercounter.cpp
	kmerpath.cpp
	probabilitycomputer.cpp
	transitionprobabilitycomputer.cpp
	uniquekmercomputer.cpp
	uniquekmers.cpp
	variant.cpp
	variantreader.cpp)

add_executable(PGGTyper main.cpp)

target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LIBRARY})
