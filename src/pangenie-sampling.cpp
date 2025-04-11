#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include "timer.hpp"
#include "commandlineparser.hpp"
#include "commands.hpp"


using namespace std;

int main(int argc, char* argv[]) {
	Timer timer;
	struct rusage rss_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping based on kmer-counting and known haplotype sequences." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	cerr << "version: v4.1.1" << endl;

	string precomputed_prefix = "";
	string readfile = "";
	string outname = "result";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	uint64_t hash_size = 3000000000;
	size_t panel_size = 0;
	double recombrate = 1.26;
	// TOD0: for testing purposes
	long double sampling_effective_N = 0.01L;
	unsigned short allele_penalty = 5;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-sampling [options] -f <index-prefix> -i <reads.fa/fq> -o <outfile-prefix>");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED");
	argument_parser.add_mandatory_argument('f', "Filename prefix of files computed by PanGenie-index (i.e. option -o used with PanGenie-index)");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_flag_argument('c', "count all read kmers instead of only those located in graph");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish");
	argument_parser.add_optional_argument('x', "0", "to which size the input panel shall be reduced.");
	argument_parser.add_optional_argument('y', "5", "Penality used for already selected alleles in sampling step.");
	argument_parser.add_optional_argument('b', "0.01", "effective population size for sampling step.");

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
	nr_jellyfish_threads = stoi(argument_parser.get_argument('j'));
	nr_core_threads = stoi(argument_parser.get_argument('t'));
	count_only_graph = !argument_parser.get_flag('c');
	panel_size = stoi(argument_parser.get_argument('x'));
	istringstream iss(argument_parser.get_argument('e'));
	iss >> hash_size;
	allele_penalty = stoi(argument_parser.get_argument('y'));
	sampling_effective_N = stof(argument_parser.get_argument('b'));


	precomputed_prefix = argument_parser.get_argument('f');

	// run sampling
	int exit_code = run_sampling(precomputed_prefix, readfile, outname, nr_jellyfish_threads, nr_core_threads, regularization, count_only_graph, hash_size, panel_size, recombrate, sampling_effective_N, allele_penalty);

	getrusage(RUSAGE_SELF, &rss_total);

	cerr << endl << "############## Summary ##############" << endl;
	cerr << "total wallclock time: \t" << timer.get_total_time() << " sec" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "#####################################" << endl;
	return exit_code;
}
