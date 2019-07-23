#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include "kmercounter.hpp"
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	clock_t clock_start = clock();
	cerr << "This is PGGTyper." << endl;
	string readfile = "";
	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	string sample_name = "sample";

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PGGTyper [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
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

	// print info
	cout << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, kmersize, sample_name);
	string segment_file = outname + "_path_segments.fasta";
	cerr << "Write path segments to file: " << segment_file << " ..." << endl;
	variant_reader.write_path_segments(segment_file);

//	// determine total genome size
//	size_t genome_kmers = variant_reader.nr_of_genomic_kmers();

	// determine chromosomes present in VCF
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);
	cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

	// TODO: only for analysis
	struct rusage r_usage0;
	getrusage(RUSAGE_SELF, &r_usage0);
	cerr << "#### Memory usage until now: " << (r_usage0.ru_maxrss / 1E6) << " GB ####" << endl;

	// determine kmer copynumbers in reads
	cerr << "Count kmers in reads ..." << endl;
	KmerCounter read_kmer_counts (readfile, kmersize);
//	cerr << "Compute kmer-coverage ..." << endl;
//	size_t kmer_coverage = read_kmer_counts.computeKmerCoverage(genome_kmers);
	size_t kmer_abundance_peak = read_kmer_counts.computeHistogram(10000, outname + "_histogram.histo");
	cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

	// count kmers in allele + reference sequence
	cerr << "Count kmers in genome ..." << endl;
	KmerCounter genomic_kmer_counts (segment_file, kmersize);

	// TODO: only for analysis
	struct rusage r_usage1;
	getrusage(RUSAGE_SELF, &r_usage1);
	cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

	// prepare output files
	variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	for (auto& chromosome : chromosomes) {
		cerr << "Processing chromosome " << chromosome << "." << endl;
		cerr << "Determine unique kmers ..." << endl;
		// determine sets of kmers unique to each variant region
		UniqueKmerComputer kmer_computer(&genomic_kmer_counts, &read_kmer_counts, &variant_reader, chromosome, kmer_abundance_peak);
		std::vector<UniqueKmers*> unique_kmers;
		kmer_computer.compute_unique_kmers(&unique_kmers);

		struct rusage r_usagei;
		getrusage(RUSAGE_SELF, &r_usagei);
		cerr << "#### Memory usage until now: " << (r_usagei.ru_maxrss / 1E6) << " GB ####" << endl;

		// get variants on this chromosome
		cerr << "Construct HMM" << endl;
		HMM hmm(&unique_kmers);

		// output the genotyping results
		cerr << "Write genotyping output ..." << endl;
		variant_reader.write_genotypes_of(chromosome, hmm.get_genotyping_result());

		// output the phasing results
		cerr << "Write phasing output ..." << endl;
		variant_reader.write_phasing_of(chromosome, hmm.get_genotyping_result());

		// destroy unique kmers
		for (size_t i = 0; i < unique_kmers.size(); ++i) {
			delete unique_kmers[i];
			unique_kmers[i] = nullptr;
		}
	}

	variant_reader.close_genotyping_outfile();
	variant_reader.close_phasing_outfile();

	cerr << endl << "###### Summary ######" << endl;
	// total time
	double cpu_time = (double)(clock() - clock_start) / CLOCKS_PER_SEC;
	cerr << "Total CPU time: " << cpu_time << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
