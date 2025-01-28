#include <iostream>
#include <sstream>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <zlib.h>
#include <cereal/archives/binary.hpp>
#include "uniquekmers.hpp"
#include "commandlineparser.hpp"



using namespace std;


struct UniqueKmersMap {
	size_t kmersize;
	std::mutex kmers_mutex;
	std::map<std::string, std::vector<std::shared_ptr<UniqueKmers>>> unique_kmers;
	std::map<std::string, double> runtimes;
	bool add_reference;

	template <class Archive>
	void save(Archive& ar) const {
		ar(kmersize, unique_kmers, runtimes, add_reference);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(kmersize, unique_kmers, runtimes, add_reference);
	}
};


struct UKAnalyzer {
	UKAnalyzer(string chromosome, shared_ptr<UniqueKmers> uk) {
		for (auto a : uk->alleles) {
			cout << chromosome << "\t" << uk->variant_pos << "\t" << a.second.kmer_path << endl;
		}
	}
};

int main (int argc, char* argv[])
{
	cerr << endl;
	cerr << "program: PanGenie - genotyping based on kmer-counting and known haplotype sequences." << endl;
	cerr << "command: Analyze-UK - Print matrix of unique kmers." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v4.0.0" << endl;

	string precomputed_uk = "";

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("Analyze-UK [options] -i <UniqueKmersMap>");
	argument_parser.add_mandatory_argument('i', "Serialized UniqueKmersMap object (.cereal).");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	precomputed_uk = argument_parser.get_argument('i');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// re-construct UniqueKmersMap from input archive
	UniqueKmersMap unique_kmers_list;
	string unique_kmers_archive = precomputed_uk;
  	ifstream os2(unique_kmers_archive, std::ios::binary);
  	cereal::BinaryInputArchive archive2( os2 );
	archive2(unique_kmers_list);

	for (auto it = unique_kmers_list.unique_kmers.begin(); it != unique_kmers_list.unique_kmers.end(); ++it) {
		for (auto uk : it->second) {
			UKAnalyzer u(it->first, uk);
		}
	}

	return 0;
}