#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <thread>
#include <algorithm>
//#include <boost/asio/thread_pool.hpp>
//#include <boost/asio/post.hpp>
//#include <boost/bind.hpp>
#include <mutex>
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

using namespace std;

/** version of the main algorithm that uses uniform transition probabilities and genotypes only based on kmer information.**/

struct Results {
	mutex result_mutex;
	map<string, vector<GenotypingResult>> result;
	map<string, double> runtimes;
};

void run_genotyping_kmers(string chromosome, KmerCounter* genomic_kmer_counts, KmerCounter* read_kmer_counts, VariantReader* variant_reader, ProbabilityTable* probs, size_t kmer_abundance_peak, bool only_genotyping, bool only_phasing, Results* results) {
	Timer timer;
	// determine sets of kmers unique to each variant region
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, variant_reader, chromosome, kmer_abundance_peak);
	std::vector<UniqueKmers*> unique_kmers;
	kmer_computer.compute_unique_kmers(&unique_kmers, probs);
	// construct HMM and run genotyping/phasing
	HMM hmm(&unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, true);
	// store the results
	{
		lock_guard<mutex> lock (results->result_mutex);
		results->result.insert(pair<string, vector<GenotypingResult>> (chromosome, move(hmm.get_genotyping_result())));
	}
	// destroy unique kmers
	for (size_t i = 0; i < unique_kmers.size(); ++i) {
		delete unique_kmers[i];
		unique_kmers[i] = nullptr;
	}
	lock_guard<mutex> lock (results->result_mutex);
	results->runtimes.insert(pair<string,double>(chromosome, timer.get_total_time()));
}

int main (int argc, char* argv[])
{
	Timer timer;
	double time_preprocessing;
	double time_kmer_counting;
	double time_writing;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie-kmers - genotyping based on kmer-counting." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	string readfile = "";
	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = false;
	bool only_phasing = false;
	long double regularization = 0.001L;;
	bool count_only_graph = true;
	bool ignore_imputed = false;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-kmers [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_flag_argument('g', "only run genotyping (Forward backward algorithm)");
	argument_parser.add_flag_argument('p', "only run phasing (Viterbi algorithm)");
	argument_parser.add_optional_argument('m', "0.001", "regularization constant for copynumber probabilities");
	argument_parser.add_flag_argument('c', "count all kmers instead of only those located in graph..");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers.");

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
	only_genotyping = argument_parser.get_flag('g');
	only_phasing = argument_parser.get_flag('p');
	regularization = stold(argument_parser.get_argument('m'));
	count_only_graph = !argument_parser.get_flag('c');
	ignore_imputed = argument_parser.get_flag('u');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, kmersize, false, sample_name);
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

	time_preprocessing = timer.get_interval_time();

	KmerCounter* read_kmer_counts = nullptr;
	// determine kmer copynumbers in reads
	if (readfile.substr(std::max(3, (int) readfile.size())-3) == std::string(".jf")) {
		cerr << "Read pre-computed read kmer counts ..." << endl;
		jellyfish::mer_dna::k(kmersize);
		read_kmer_counts = new JellyfishReader(readfile, kmersize);
	} else {
		cerr << "Count kmers in reads ..." << endl;
		if (count_only_graph) {
			read_kmer_counts = new JellyfishCounter(readfile, segment_file, kmersize, nr_jellyfish_threads);
		} else {
			read_kmer_counts = new JellyfishCounter(readfile, kmersize, nr_jellyfish_threads);
		}
	}

	size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
	cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

	// count kmers in allele + reference sequence
	cerr << "Count kmers in genome ..." << endl;

	JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads);

	// TODO: only for analysis
	struct rusage r_usage1;
	getrusage(RUSAGE_SELF, &r_usage1);
	cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

	// prepare output files
	if (! only_phasing) variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	if (! only_genotyping) variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	time_kmer_counting = timer.get_interval_time();

	cerr << "Construct HMM and run core algorithm ..." << endl;
	ProbabilityTable probabilities (kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);

	// determine max number of available threads (at most one thread per chromosome possible)
	size_t available_threads = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
	if (nr_core_threads > available_threads) {
		cerr << "Warning: set nr_core_threads to " << available_threads << "." << endl;
		nr_core_threads = available_threads;
 	}
	Results results;
	{
		// create thread pool
		ThreadPool threadPool (nr_core_threads);
		for (auto chromosome : chromosomes) {
			KmerCounter* genomic = &genomic_kmer_counts;
			VariantReader* variants = &variant_reader;
			Results* r = &results;
			ProbabilityTable* probs = &probabilities;
			function<void()> f_genotyping = bind(run_genotyping_kmers, chromosome, genomic, read_kmer_counts, variants, probs, kmer_abundance_peak, only_genotyping, only_phasing, r);
			threadPool.submit(f_genotyping);
		}
	}

/**
	// create thread pool
	boost::asio::thread_pool threadPool(nr_core_threads);
	for (auto chromosome : chromosomes) {
		boost::asio::post(threadPool, boost::bind(run_genotyping_kmers, chromosome, &genomic_kmer_counts, read_kmer_counts, &variant_reader, kmer_abundance_peak, only_genotyping, only_phasing, &results));
	} 
	threadPool.join();
**/

	timer.get_interval_time();

	// output VCF
	cerr << "Write results to VCF ..." << endl;
	assert (results.result.size() == chromosomes.size());
	// write VCF
	for (auto it = results.result.begin(); it != results.result.end(); ++it) {
		if (!only_phasing) {
			// output genotyping results
			variant_reader.write_genotypes_of(it->first, it->second, ignore_imputed);
		}
		if (!only_genotyping) {
			// output phasing results
			variant_reader.write_phasing_of(it->first, it->second, ignore_imputed);
		}
	}

	if (! only_phasing) variant_reader.close_genotyping_outfile();
	if (! only_genotyping) variant_reader.close_phasing_outfile();

	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers:\t" << time_kmer_counting << " sec" << endl;
	// output per chromosome time
	double time_hmm = time_writing;
	for (auto chromosome : chromosomes) {
		double time_chrom = results.runtimes.at(chromosome);
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << time_chrom << endl;
		time_hmm += time_chrom;
	}
	cerr << "total running time:\t" << time_preprocessing + time_kmer_counting + time_hmm << " sec"<< endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	delete read_kmer_counts;
	return 0;
}

