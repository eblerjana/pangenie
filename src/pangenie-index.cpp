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
#include "uniquekmercomputer.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "threadpool.hpp"



using namespace std;

bool ends_with (string const &filename, string const &ending) {
    if (filename.length() >= ending.length()) {
        return (0 == filename.compare (filename.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void check_input_file(string &filename) {
	// check if file exists and can be opened
	ifstream file(filename);
	if (!file.good()) {
		stringstream ss;
		ss << "File " << filename << " cannot be opened." << endl;
		throw runtime_error(ss.str());
	}
	// make sure file is not compressed
	if (ends_with(filename, ".gz")) {
		stringstream ss;
		ss << "File " << filename << " seems to be gzip-compressed. PanGenie requires an uncompressed file." << endl;
		throw runtime_error(ss.str());
	}
}

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


void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, VariantReader* variant_reader, UniqueKmersMap* unique_kmers_map, string outname) {
	Timer timer;
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, variant_reader, chromosome);
	std::vector<shared_ptr<UniqueKmers>> unique_kmers;
	string filename = outname + "_" + chromosome + "_kmers.tsv.gz";
	kmer_computer.compute_unique_kmers(&unique_kmers, filename, true);
	// store the results
	{
		lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
		unique_kmers_map->unique_kmers.insert(pair<string, vector<shared_ptr<UniqueKmers>>> (chromosome, move(unique_kmers)));
	}
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes.insert(pair<string, double>(chromosome, timer.get_total_time()));
}


int main (int argc, char* argv[])
{
	Timer timer;
	double time_preprocessing;
	double time_graph_counting;
	double time_unique_kmers;
	double time_writing;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping based on kmer-counting and known haplotype sequences." << endl;
	cerr << "command: PanGenie-index - construct graph and determine unique kmers." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v3.0.0" << endl;

	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	size_t nr_jellyfish_threads = 1;
	bool add_reference = true;
	uint64_t hash_size = 3000000000;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-index [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders.");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish.");
//	argument_parser.add_flag_argument('d', "do not add reference as additional path.");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	reffile = argument_parser.get_argument('r');
	vcffile = argument_parser.get_argument('v');
	kmersize = stoi(argument_parser.get_argument('k'));
	outname = argument_parser.get_argument('o');
	nr_jellyfish_threads = stoi(argument_parser.get_argument('t'));
	istringstream iss(argument_parser.get_argument('e'));
	iss >> hash_size;
	add_reference = !argument_parser.get_flag('d');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// check if input files exist and are uncompressed
	check_input_file(reffile);
	check_input_file(vcffile);

	UniqueKmersMap unique_kmers_list;
	vector<string> chromosomes;
	string segment_file = outname + "_path_segments.fasta";
	size_t available_threads_uk;
	size_t nr_cores_uk;
	unsigned short nr_paths;

	/**
	*  1) Indexing step. Read variant information and determine unique kmers.
	*/

	{
		/** 
		*   Step 1: read variants, merge variants that are closer than kmersize apart,
		*   and write allele sequences and unitigs inbetween to a file.
		**/ 
		cerr << "Determine allele sequences ..." << endl;
		VariantReader variant_reader (vcffile, reffile, kmersize, add_reference);

		cerr << "Write path segments to file: " << segment_file << " ..." << endl;
		variant_reader.write_path_segments(segment_file);

		// determine chromosomes present in VCF
		variant_reader.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		// print RSS up to now
		struct rusage r_usage00;
		getrusage(RUSAGE_SELF, &r_usage00);
		cerr << "#### Max RSS after determing allele sequences: " << (r_usage00.ru_maxrss / 1E6) << " GB ####" << endl;

		time_preprocessing = timer.get_interval_time();
		nr_paths = variant_reader.nr_of_paths();

	
		/**
		*  Step 2: count graph k-mers. Needed to determine unique k-mers in subsequent steps.
		**/ 
		cerr << "Count kmers in graph ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);

		// print RSS up to now
		struct rusage r_usage1;
		getrusage(RUSAGE_SELF, &r_usage1);
		cerr << "#### Max RSS after counting graph kmers: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

		time_graph_counting = timer.get_interval_time();


		/**
		* Step 3: determine unique k-mers for each variant bubble and prepare datastructure storing
		* that information. It will later be used for the genotyping step. 
		* This step will delete information from the VariantReader that will no longer be used.
		**/

		cerr << "Determine unique kmers ..." << endl;
		available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		nr_cores_uk = min(nr_jellyfish_threads, available_threads_uk);
		if (nr_cores_uk < nr_jellyfish_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}


		{
			// create thread pool with at most nr_chromosome threads
			ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				VariantReader* variants = &variant_reader;
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, variants, result, outname);
				threadPool.submit(f_unique_kmers);
			}
		}

		// print RSS up to now
		struct rusage r_usage2;
		getrusage(RUSAGE_SELF, &r_usage2);
		cerr << "#### Max RSS after determing unique kmers: " << (r_usage2.ru_maxrss / 1E6) << " GB ####" << endl;

		// determine the total runtime needed to compute unique kmers
		time_unique_kmers = 0.0;
		for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
			time_unique_kmers += it->second;
		}
		timer.get_interval_time();	
	}

	// serialization of UniqueKmersMap object
	cerr << "Storing unique kmer information ..." << endl;
  	ofstream os(outname + "_UniqueKmersMap.cereal", std::ios::binary);
  	cereal::BinaryOutputArchive archive( os );
	archive(unique_kmers_list);
	time_writing = timer.get_interval_time();

	// print RSS up to now
	struct rusage r_usage3;
	getrusage(RUSAGE_SELF, &r_usage3);
	cerr << "#### Max RSS after Indexing is done and objects are deleted: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;

	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting graph kmers (wallclock):\t" << time_graph_counting << " sec" << endl;
	cerr << "time spent determining unique kmers:\t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent writing output file:\t" << time_writing << " sec" << endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;


	return 0;
}
