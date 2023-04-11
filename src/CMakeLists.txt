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
	stepwiseuniquekmercomputer.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	uniquekmers.cpp
	variant.cpp
	variantreader.cpp)

add_executable(PanGenie pangenie.cpp)
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})

add_executable(PanGenie-serialize pangenie-serialize.cpp)
target_link_libraries(PanGenie-serialize PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-serialize PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})

add_executable(PanGenie-index pangenie-index.cpp)
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})

add_executable(PanGenie-genotype pangenie-genotype.cpp)
target_link_libraries(PanGenie-genotype PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-genotype PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})

add_executable(Analyze analyze.cpp)
target_link_libraries(Analyze PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(Analyze PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})