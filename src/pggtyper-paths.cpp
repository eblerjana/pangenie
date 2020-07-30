#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <thread>
#include <algorithm>
//#include <boost/asio/thread_pool.hpp>
//#include <boost/asio/post.hpp>
//#include <boost/bind.hpp>
#include <mutex>
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "threadpool.hpp"
#include <cassert>

using namespace std;

// TODO
/** version of the main algorithm that uses only path information for genotyping and ignores kmers **/

struct Results {
	mutex result_mutex;
	map<string, vector<GenotypingResult>> result;
	map<string, double> runtimes;
};

void run_genotyping_paths(string chromosome, KmerCounter* genomic_kmer_counts, KmerCounter* read_kmer_counts, VariantReader* variant_reader, ProbabilityTable* probs, size_t kmer_abundance_peak, bool only_genotyping, bool only_phasing, long double effective_N, Results* results) {
	Timer timer;
	// determine sets of kmers unique to each variant region
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, variant_reader, chromosome, kmer_abundance_peak);
	std::vector<UniqueKmers*> unique_kmers;
	kmer_computer.compute_empty(&unique_kmers);
	// construct HMM and run genotyping/phasing
	HMM hmm(&unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, false, effective_N);
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
	double time_writing;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie-paths - genotyping and phasing based on known haplotype paths." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	string reffile = "";
	string vcffile = "";
	string outname = "result";
	string sample_name = "sample";
	size_t nr_core_threads = 1;
	bool only_genotyping = false;
	bool only_phasing = false;
	double effective_N = 0.00001L;
	bool add_reference = true;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie-paths [options] -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_optional_argument('n', "0.00001", "effective population size");
	argument_parser.add_flag_argument('g', "only run genotyping (Forward backward algorithm)");
	argument_parser.add_flag_argument('p', "only run phasing (Viterbi algorithm)");
	argument_parser.add_flag_argument('d', "do not add reference as additional path.");

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
	nr_core_threads = stoi(argument_parser.get_argument('t'));
	only_genotyping = argument_parser.get_flag('g');
	only_phasing = argument_parser.get_flag('p');
	effective_N = stold(argument_parser.get_argument('n'));
	add_reference = !argument_parser.get_flag('d');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, 31, add_reference, sample_name);

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

	cerr << "Construct HMM and run core algorithm ..." << endl;
	ProbabilityTable probabilities;

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
			KmerCounter* genomic = nullptr;
			VariantReader* variants = &variant_reader;
			Results* r = &results;
			ProbabilityTable* probs = &probabilities;
			function<void()> f_genotyping = bind(run_genotyping_paths, chromosome, genomic, nullptr, variants, probs, 28, only_genotyping, only_phasing, effective_N, r);
			threadPool.submit(f_genotyping);
		}
	}

/**	boost::asio::thread_pool threadPool(nr_core_threads);
	for (auto chromosome : chromosomes) {
		boost::asio::post(threadPool, boost::bind(run_genotyping_paths, chromosome, nullptr, nullptr, &variant_reader, 28, only_genotyping, only_phasing, effective_N, &results));
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
			variant_reader.write_genotypes_of(it->first, it->second);
		}
		if (!only_genotyping) {
			// output phasing results
			variant_reader.write_phasing_of(it->first, it->second);
		}
	}

	if (! only_phasing) variant_reader.close_genotyping_outfile();
	if (! only_genotyping) variant_reader.close_phasing_outfile();

	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	// output per chromosome time
	double time_hmm = time_writing;
	for (auto chromosome : chromosomes) {
		double time_chrom = results.runtimes.at(chromosome);
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << time_chrom << endl;
		time_hmm += time_chrom;
	}
	cerr << "total running time:\t" << time_preprocessing + time_hmm << " sec"<< endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}

