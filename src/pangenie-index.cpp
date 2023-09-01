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
	cerr << "version: v3.0.1" << endl;

	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	size_t nr_jellyfish_threads = 1;
	bool add_reference = true;
	uint64_t hash_size = 3000000000;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-index [options] -r <reference.fa> -v <variants.vcf> -o <index-prefix>");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED");
	argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED");
	argument_parser.add_mandatory_argument('o', "prefix of the output files. NOTE: the given path must not include non-existent folders");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish");
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

	// run preprocessing
	int exit_code = run_index_command(reffile, vcffile, kmersize, outname, nr_jellyfish_threads, add_reference, hash_size);
	getrusage(RUSAGE_SELF, &rss_total);


	cerr << endl << "############## Summary ##############" << endl;
	cerr << "total wallclock time: \t" << timer.get_total_time() << " sec" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "#####################################" << endl;

	return exit_code;
}