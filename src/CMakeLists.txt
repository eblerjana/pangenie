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

add_executable(PanGenie pggtyper.cpp)
add_executable(PanGenie-kmers pggtyper-kmers.cpp)
add_executable(PanGenie-paths pggtyper-paths.cpp)
add_executable(PanGenie-graph pggtyper-graph.cpp)

#target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LDFLAGS_OTHER})
#target_link_libraries(PGGTyper PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LIBRARIES})
#target_link_libraries(PGGTyper-kmers PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PanGenie-kmers PanGenieLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie-kmers PanGenieLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LIBRARIES})
#target_link_libraries(PGGTyper-paths PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PanGenie-paths PanGenieLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie-paths PanGenieLib ${JELLYFISH_LIBRARIES})

#target_link_libraries(PGGTyper-graph PGGTyperLib ${JELLYFISH_LDFLAGS_OTHER} ${Boost_LIBRARIES})
#target_link_libraries(PGGTyper-graph PGGTyperLib ${JELLYFISH_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(PanGenie-graph PanGenieLib ${JELLYFISH_LDFLAGS_OTHER})
target_link_libraries(PanGenie-graph PanGenieLib ${JELLYFISH_LIBRARIES})
