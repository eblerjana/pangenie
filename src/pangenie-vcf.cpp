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
	cerr << "version: v4.2.1" << endl;

	string precomputed_prefix = "";
	string outname = "result";
	string sample_name = "sample";
	string results_name = "";

	bool only_genotyping = true;
	bool only_phasing = false;
	bool ignore_imputed = false;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-vcf [options] -f <index-prefix> -z <outname_genotyping.cereal> -o <outfile-prefix>");
	argument_parser.add_mandatory_argument('z', "serialized genotyping results (produced by PanGenie run with parameter -w)");
	argument_parser.add_mandatory_argument('f', "filename prefix of the index files (i.e. option -o used with PanGenie-index)");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_flag_argument('g', "run genotyping (only set if used with PanGenie genotyping)");
	argument_parser.add_flag_argument('p', "run phasing (only set if used with PanGenie genotyping)");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers");


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

	outname = argument_parser.get_argument('o');
	sample_name = argument_parser.get_argument('s');
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

	ignore_imputed = argument_parser.get_flag('u');
	precomputed_prefix = argument_parser.get_argument('f');
	results_name = argument_parser.get_argument('z');

	int exit_code = run_vcf_command(precomputed_prefix, results_name, outname, sample_name, only_genotyping, only_phasing, ignore_imputed);
	getrusage(RUSAGE_SELF, &rss_total);
	cerr << endl << "############## Summary ##############" << endl;
	cerr << "total wallclock time: \t" << timer.get_total_time() << " sec" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "#####################################" << endl;
	return exit_code;
}
