add_library(PGGTyperLib SHARED 
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	commandlineparser.cpp
	columnindexer.cpp
	copynumberassignment.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	histogram.cpp
	hmm.cpp
	jellyfishcounter.cpp
	jellyfishreader.cpp
	kmerpath.cpp
	probabilitycomputer.cpp
	sequenceutils.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	uniquekmers.cpp
	variant.cpp
	variantreader.cpp)

add_executable(PGGTyper pggtyper.cpp)
add_executable(PGGTyper-kmers pggtyper-kmers.cpp)
add_executable(PGGTyper-paths pggtyper-paths.cpp)

#target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LDFLAGS_OTHER})
#target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LIBRARIES})
#target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LIBRARIES})
#target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LIBRARIES})
