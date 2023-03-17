#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include "kmercounter.hpp"
#include "jellyfishreader.hpp"
#include "jellyfishcounter.hpp"
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "threadpool.hpp"
#include "pathsampler.hpp"

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
	map<string, vector<UniqueKmers*>> unique_kmers;
	map<string, double> runtimes;
};

int main (int argc, char* argv[])
{
	Timer timer;
	double time_preprocessing;
	double time_graph_counting;
	double time_unique_kmers;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping and phasing based on kmer-counting and known haplotype sequences." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v3.0.0" << endl;
	string readfile = "";
	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.001L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	bool add_reference = true;
	size_t sampling_size = 0;
	uint64_t hash_size = 3000000000;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders.");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_flag_argument('g', "run genotyping (Forward backward algorithm, default behaviour).");
	argument_parser.add_flag_argument('p', "run phasing (Viterbi algorithm). Experimental feature.");
	argument_parser.add_flag_argument('c', "count all read kmers instead of only those located in graph.");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers.");
	argument_parser.add_flag_argument('d', "do not add reference as additional path.");
	argument_parser.add_optional_argument('a', "0", "sample subsets of paths of this size.");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish.");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}
	readfile = argument_parser.get_argument('i');
	reffile = argument_parser.get_argument('r');
	vcffile = argument_parser.get_argument('v');
	kmersize = stoi(argument_parser.get_argument('k'));
	outname = argument_parser.get_argument('o');
	sample_name = argument_parser.get_argument('s');
	nr_jellyfish_threads = stoi(argument_parser.get_argument('j'));
	nr_core_threads = stoi(argument_parser.get_argument('t'));

	bool genotyping_flag = argument_parser.get_flag('g');
	bool phasing_flag = argument_parser.get_flag('p');
	
	if (genotyping_flag && phasing_flag) {
		only_genotyping = false;
		only_genotyping = false;
	}
	if (!genotyping_flag && phasing_flag) {
		only_genotyping = false;
		only_phasing = true;
	}

	count_only_graph = !argument_parser.get_flag('c');
	ignore_imputed = argument_parser.get_flag('u');
	add_reference = !argument_parser.get_flag('d');
	sampling_size = stoi(argument_parser.get_argument('a'));
	istringstream iss(argument_parser.get_argument('e'));
	iss >> hash_size;

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// check if input files exist and are uncompressed
	check_input_file(reffile);
	check_input_file(vcffile);
	check_input_file(readfile);

	UniqueKmersMap unique_kmers_list;

	/**
	*  1) Indexing step. Read variant information and determine unique kmers.
	*/

	{
		/** 
		*   Step 1: read variants, merge variants that are closer than kmersize apart,
		*   and write allele sequences and unitigs inbetween to a file.
		*   TODO: store Variant objects as pointers and implement function to delete them (in order to save memory)
		**/ 
		cerr << "Determine allele sequences ..." << endl;
		VariantReader variant_reader (vcffile, reffile, kmersize, add_reference, sample_name);

		string segment_file = outname + "_path_segments.fasta";
		cerr << "Write path segments to file: " << segment_file << " ..." << endl;
		variant_reader.write_path_segments(segment_file);

		// determine chromosomes present in VCF
		vector<string> chromosomes;
		variant_reader.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		// print RSS up to now
		struct rusage r_usage00;
		getrusage(RUSAGE_SELF, &r_usage00);
		cerr << "#### Max RSS after determing allele sequences: " << (r_usage00.ru_maxrss / 1E6) << " GB ####" << endl;

		time_preprocessing = timer.get_interval_time();


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
		* TODO: store actual k-mer sequences and leave an empty structure for counts to be filled later once read-kmer counts are available
		* TODO: in old version, kmers were selected taking kmer_coverage into account. This is no longer possible, because the step
		*       counting kmers in the reads happens later.
		* TODO: write an output file which lists unique kmers per variant, as well as left + right overhangs (needed in subsequent steps + usefull information anyways)
		*       the UniqueKmers object therefore does not need to store any unique kmers or overhang information
		**/

		cerr << "Determine unique kmers ..." << endl;
		size_t available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		size_t nr_cores_uk = min(nr_core_threads, available_threads_uk);
		if (nr_cores_uk < nr_core_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}


		{
			// create thread pool with at most nr_chromosome threads
			ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				VariantReader* variants = &variant_reader;
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, variants, result);
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

	// print RSS up to now
	struct rusage r_usage3;
	getrusage(RUSAGE_SELF, &r_usage3);
	cerr << "#### Max RSS after Indexing is done and objects are deleted: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;

	return 0;
}
