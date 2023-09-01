#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include "timer.hpp"
#include "commands.hpp"
#include "commandlineparser.hpp"


using namespace std;

int main(int argc, char* argv[]) {
	Timer timer;
	struct rusage rss_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping based on kmer-counting and known haplotype sequences." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v3.0.1" << endl;

	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;

	string precomputed_prefix = "";
	string readfile = "";
	string outname = "result";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 0;
	uint64_t hash_size = 3000000000;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie [options] -f <index-prefix> -i <reads.fa/fq> -o <outfile-prefix>\nPanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -o <outfile-prefix>");
	argument_parser.add_optional_argument('r', "", "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED");
	argument_parser.add_optional_argument('v', "", "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED");
	argument_parser.add_optional_argument('f', "", "Filename prefix of files computed by PanGenie-index (i.e. option -o used with PanGenie-index)");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_flag_argument('g', "run genotyping (Forward backward algorithm, default behaviour)");
	argument_parser.add_flag_argument('p', "run phasing (Viterbi algorithm). Experimental feature");
	argument_parser.add_flag_argument('c', "count all read kmers instead of only those located in graph");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers");
	argument_parser.add_optional_argument('a', "0", "sample subsets of paths of this size");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish");

	argument_parser.exactly_one('f', 'v');
	argument_parser.exactly_one('f', 'r');
	argument_parser.not_both('f', 'k');

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	readfile = argument_parser.get_argument('i');
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
	sampling_size = stoi(argument_parser.get_argument('a'));
	istringstream iss(argument_parser.get_argument('e'));
	iss >> hash_size;

	if (argument_parser.exists('f')) {
		precomputed_prefix = argument_parser.get_argument('f');

		// run genotyping
		int exit_code = run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size);

		getrusage(RUSAGE_SELF, &rss_total);


		cerr << endl << "############## Summary ##############" << endl;
		cerr << "total wallclock time: \t" << timer.get_total_time() << " sec" << endl;
		cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
		cerr << "#####################################" << endl;
		return exit_code;
	} else {

		reffile = argument_parser.get_argument('r');
		vcffile = argument_parser.get_argument('v');
		kmersize = stoi(argument_parser.get_argument('k'));
		bool add_reference = true;

		cerr << endl << "NOTE: by running PanGenie-index first to pre-process data, you can reduce memory usage and speed up PanGenie. This is helpful especially when genotyping the same variants across multiple samples." << endl << endl;

		int exit_code = run_single_command(outname, readfile, reffile, vcffile, kmersize, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, add_reference, sampling_size, hash_size);

		getrusage(RUSAGE_SELF, &rss_total);

		cerr << endl << "############## Summary ##############" << endl;
		cerr << "total wallclock time: \t" << timer.get_total_time() << " sec" << endl;
		cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
		cerr << "#####################################" << endl;
		return exit_code;
	}
}