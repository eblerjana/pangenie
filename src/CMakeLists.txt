add_library(PanGenieLib SHARED 
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	commandlineparser.cpp
	commands.cpp
	columnindexer.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	graphbuilder.cpp
	graph.cpp
	haplotypesampler.cpp
	histogram.cpp
	hmm.cpp
	jellyfishcounter.cpp
	jellyfishreader.cpp
	kmerpath.cpp
	kmerpath16.cpp
	kmerparser.cpp
	pathsampler.cpp
	probabilitycomputer.cpp
	probabilitytable.cpp
	sampledpanel.cpp
	samplingemissions.cpp
	samplingtransitions.cpp
	sequenceutils.cpp
	stepwiseuniquekmercomputer.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	multiallelicuniquekmers.cpp
	biallelicuniquekmers.cpp
	variant.cpp
	variantreader.cpp)


add_executable(PanGenie-index pangenie-index.cpp)
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie pangenie-genotype.cpp)
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie-vcf pangenie-vcf.cpp)
target_link_libraries(PanGenie-vcf PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-vcf PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie-sampling pangenie-sampling.cpp)
target_link_libraries(PanGenie-sampling PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-sampling PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(Analyze-UK analyze-uk.cpp)
target_link_libraries(Analyze-UK PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(Analyze-UK PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})
