#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"

using namespace std;

// TODO
/** version of the main algorithm that uses only path information for genotyping and ignores kmers **/

int main (int argc, char* argv[])
{
	Timer timer;
	double time_preprocessing;
	double time_total;

	cerr << endl;
	cerr << "program: PGGTyper-paths - genotyping and phasing based on known haplotype paths." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	string reffile = "";
	string vcffile = "";
	string outname = "result";
	string sample_name = "sample";
	bool only_genotyping = false;
	bool only_phasing = false;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PGGTyper-paths [options] -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_flag_argument('g', "only run genotyping (Forward backward algorithm).");
	argument_parser.add_flag_argument('p', "only run phasing (Viterbi algorithm).");
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
	outname = argument_parser.get_argument('o');
	sample_name = argument_parser.get_argument('s');
	only_genotyping = argument_parser.get_flag('g');
	only_phasing = argument_parser.get_flag('p');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, 31, sample_name);
	string segment_file = outname + "_path_segments.fasta";
	cerr << "Write path segments to file: " << segment_file << " ..." << endl;
	variant_reader.write_path_segments(segment_file);

	// determine chromosomes present in VCF
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);
	cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

	// TODO: only for analysis
	struct rusage r_usage0;
	getrusage(RUSAGE_SELF, &r_usage0);
	cerr << "#### Memory usage until now: " << (r_usage0.ru_maxrss / 1E6) << " GB ####" << endl;

	// prepare output files
	if (! only_phasing) variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	if (! only_genotyping) variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	time_preprocessing = timer.get_interval_time();

	map<string,double> time_hmm;
	for (auto& chromosome : chromosomes) {
		Timer timer_chrom;
		cerr << "Processing chromosome " << chromosome << "." << endl;
		UniqueKmerComputer kmer_computer(nullptr, nullptr, &variant_reader, chromosome, 28);
		std::vector<UniqueKmers*> unique_kmers;
		kmer_computer.compute_empty(&unique_kmers);

		struct rusage r_usagei;
		getrusage(RUSAGE_SELF, &r_usagei);
		cerr << "#### Memory usage until now: " << (r_usagei.ru_maxrss / 1E6) << " GB ####" << endl;

		// get variants on this chromosome
		cerr << "Construct HMM" << endl;
		HMM hmm(&unique_kmers, !only_phasing, !only_genotyping);

		if (! only_phasing) {
			// output the genotyping results
			cerr << "Write genotyping output ..." << endl;
			variant_reader.write_genotypes_of(chromosome, hmm.get_genotyping_result());
		}

		if (! only_genotyping) {
			// output the phasing results
			cerr << "Write phasing output ..." << endl;
			variant_reader.write_phasing_of(chromosome, hmm.get_genotyping_result());
		}

		// destroy unique kmers
		for (size_t i = 0; i < unique_kmers.size(); ++i) {
			delete unique_kmers[i];
			unique_kmers[i] = nullptr;
		}
		double time_chrom = timer_chrom.get_total_time();
		time_hmm[chromosome] = time_chrom;
	}

	if (! only_phasing) variant_reader.close_genotyping_outfile();
	if (! only_genotyping) variant_reader.close_phasing_outfile();

	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	for (auto chromosome : chromosomes) {
		cerr << "time spent genotyping chromosome " << chromosome << ": " << time_hmm.at(chromosome) << endl; 
	}
	cerr << "total wallclock time: " << time_total  << " sec" << endl;
	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
