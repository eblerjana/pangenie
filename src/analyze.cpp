#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <zlib.h>
#include <cereal/archives/binary.hpp>
#include "kmercounter.hpp"
#include "jellyfishreader.hpp"
#include "jellyfishcounter.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "stepwiseuniquekmercomputer.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "hmm.hpp"
#include "threadpool.hpp"



using namespace std;

struct UniqueKmersMap {
	mutex kmers_mutex;
	map<string, vector<shared_ptr<UniqueKmers>>> unique_kmers;
	map<string, double> runtimes;

	template <class Archive>
	void save(Archive& ar) const {
		ar(unique_kmers, runtimes);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(unique_kmers, runtimes);
	}
};

struct Results {
	mutex result_mutex;
	map<string, shared_ptr<HMM>> result;
	map<string, double> runtimes;

	template <class Archive>
	void save(Archive& ar) const {
		ar(result, runtimes);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(result, runtimes);
	}
};


int main (int argc, char* argv[])
{
	Timer timer;

	struct rusage rss_variant_reader;
	struct rusage rss_unique_kmers_map;
	struct rusage rss_results;

	cerr << endl;
	cerr << "program: PanGenie - genotyping based on kmer-counting and known haplotype sequences." << endl;
	cerr << "command: PanGenie-index - construct graph and determine unique kmers." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v3.0.0" << endl;

	string precomputed_prefix = "";

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("Analyze [options] -i <prefix>");
	argument_parser.add_mandatory_argument('i', "Prefix of the .cereal Files produced by PanGenie-index and PanGenie-genotype.");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	precomputed_prefix = argument_parser.get_argument('i');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// re-construct VariantReader from input archive
	VariantReader variant_reader;
	string variant_reader_archive = precomputed_prefix + "_VariantReader.cereal";
	cerr << "Reading precomputed VariantReader from " << variant_reader_archive << " ..." << endl; 
  	ifstream os1(variant_reader_archive, std::ios::binary);
  	cereal::BinaryInputArchive archive1( os1 );
	archive1(variant_reader);
	getrusage(RUSAGE_SELF, &rss_variant_reader);


	// re-construct UniqueKmersMap from input archive
	UniqueKmersMap unique_kmers_list;
	string unique_kmers_archive = precomputed_prefix + "_UniqueKmersMap.cereal";
	cerr << "Reading precomputed UniqueKmersMap from " << unique_kmers_archive << " ..." << endl; 
  	ifstream os2(unique_kmers_archive, std::ios::binary);
  	cereal::BinaryInputArchive archive2( os2 );
	archive2(unique_kmers_list);
	getrusage(RUSAGE_SELF, &rss_unique_kmers_map);

	// re-construct Results from input archive
	Results results;
	string results_archive = precomputed_prefix + "_Results.cereal";
	cerr << "Reading precomputed Results from " << results_archive << " ..." << endl; 
  	ifstream os3(results_archive, std::ios::binary);
  	cereal::BinaryInputArchive archive3( os3 );
	archive3(results);
	getrusage(RUSAGE_SELF, &rss_results);


	// print max RSS after each step
	cerr << "Max RSS after reading VariantReader from disk: \t" << (rss_variant_reader.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after reading UniqueKmersMap from disk: \t" << (rss_unique_kmers_map.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after reading Results from disk: \t" << (rss_results.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
