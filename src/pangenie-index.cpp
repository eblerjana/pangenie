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
#include "graphbuilder.hpp"
#include "stepwiseuniquekmercomputer.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "threadpool.hpp"
#include "graph.hpp"



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


void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, shared_ptr<Graph> graph, UniqueKmersMap* unique_kmers_map, string outname) {
	Timer timer;
	StepwiseUniqueKmerComputer kmer_computer(genomic_kmer_counts, graph);
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
	double time_preprocessing = 0.0;
	double time_kmer_counting = 0.0;
	double time_serialize_graph = 0.0;
	double time_unique_kmers = 0.0;
	double time_serialize = 0.0;
	double time_total = 0.0;

	struct rusage rss_preprocessing;
	struct rusage rss_kmer_counting;
	struct rusage rss_serialize_graph;
	struct rusage rss_unique_kmers;
	struct rusage rss_total;

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
	argument_parser.add_command("PanGenie-index [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -o <outfile-prefix>");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('o', "prefix of the output files. NOTE: the given path must not include non-existent folders.");
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
//	add_reference = !argument_parser.get_flag('d');

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

	/**
	*  1) Indexing step. Read variant information and determine unique kmers.
	*/

	{
		/** 
		*   Step 1: read variants, merge variants that are closer than kmersize apart,
		*   and write allele sequences and unitigs inbetween to a file.
		**/ 
		cerr << "Determine allele sequences ..." << endl;
		map<string, shared_ptr<Graph>> graph;
		GraphBuilder graph_builder (vcffile, reffile, graph, segment_file, kmersize, add_reference);

		// determine chromosomes present in VCF
		graph_builder.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		getrusage(RUSAGE_SELF, &rss_preprocessing);
		time_preprocessing = timer.get_interval_time();

		/**
		*  Step 2: count graph k-mers. Needed to determine unique k-mers in subsequent steps.
		**/ 
		cerr << "Count kmers in graph ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);


		getrusage(RUSAGE_SELF, &rss_kmer_counting);
		time_kmer_counting = timer.get_interval_time();

		cerr << "Serialize Graph objects ..." << endl;
		for (auto chromosome : chromosomes) {
  			ofstream os(outname + "_" + chromosome + "_Graph.cereal", std::ios::binary);
  			cereal::BinaryOutputArchive archive( os );
			archive(*graph.at(chromosome));;
		}

		getrusage(RUSAGE_SELF, &rss_serialize_graph);
		time_serialize_graph = timer.get_interval_time();


		/**
		* Step 3: determine unique k-mers for each variant bubble and prepare datastructure storing
		* that information. It will later be used for the genotyping step. 
		* This step will delete information from the GraphBuilder that will no longer be used.
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
				shared_ptr<Graph> graph_segment = graph.at(chromosome);
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, graph_segment, result, outname);
				threadPool.submit(f_unique_kmers);
			}
		}

		// determine the total runtime needed to compute unique kmers
		time_unique_kmers = 0.0;
		for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
			time_unique_kmers += it->second;
		}

		getrusage(RUSAGE_SELF, &rss_unique_kmers);
		timer.get_interval_time();	
	}

	// serialization of UniqueKmersMap object
	cerr << "Storing unique kmer information ..." << endl;
	{
  		ofstream os(outname + "_UniqueKmersMap.cereal", std::ios::binary);
  		cereal::BinaryOutputArchive archive( os );
		archive(unique_kmers_list);
	}

	getrusage(RUSAGE_SELF, &rss_total);
	time_serialize = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers in genome (wallclock): \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent writing Graph objects to disk: \t" << time_serialize_graph << " sec" << endl;
	cerr << "time spent determining unique kmers: \t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent writing UniqueKmersMap to disk: \t" << time_serialize << " sec" << endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	cerr << endl;
	cerr << "Max RSS after reading input files: \t" << (rss_preprocessing.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers in genome: \t" << (rss_kmer_counting.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after determining unique kmers: \t" << (rss_unique_kmers.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
	return 0;
}