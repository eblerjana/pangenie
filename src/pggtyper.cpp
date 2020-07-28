#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <mutex>
#include <thread>
#include <algorithm>
//#include <boost/asio/thread_pool.hpp>
//#include <boost/asio/post.hpp>
//#include <boost/bind.hpp>
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

struct UniqueKmersMap {
	mutex kmers_mutex;
	map<string, vector<UniqueKmers*>> unique_kmers;
	map<string, double> runtimes;
};

struct Results {
	mutex result_mutex;
	map<string, vector<GenotypingResult>> result;
	map<string, double> runtimes;
};

// bind(prepare_unique_kmers, chromosome, genomic, read_kmer_counts, variants, kmer_abundance_peak, result);
void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, KmerCounter* read_kmer_counts, VariantReader* variant_reader, size_t kmer_abundance_peak, long double regularization, UniqueKmersMap* unique_kmers_map) {
	Timer timer;
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, variant_reader, chromosome, kmer_abundance_peak);
	std::vector<UniqueKmers*> unique_kmers;
	kmer_computer.compute_unique_kmers(&unique_kmers, regularization);
	// store the results
	{
		lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
		unique_kmers_map->unique_kmers.insert(pair<string, vector<UniqueKmers*>> (chromosome, move(unique_kmers)));
	}
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes.insert(pair<string, double>(chromosome, timer.get_total_time()));
}

void run_genotyping(string chromosome, vector<UniqueKmers*>* unique_kmers, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, vector<size_t>* only_paths, Results* results) {
	Timer timer;
	// construct HMM and run genotyping/phasing
	HMM hmm(unique_kmers, !only_phasing, !only_genotyping, 1.26, false, effective_N, only_paths, false);
	// store the results
	{
		lock_guard<mutex> lock_result (results->result_mutex);
		// combine the new results to the already existing ones (if present)
		if (results->result.find(chromosome) == results->result.end()) {
			results->result.insert(pair<string, vector<GenotypingResult>> (chromosome, move(hmm.get_genotyping_result())));
		} else {
			// combine newly computed likelihoods with already exisiting ones
			size_t index = 0;
			for (auto likelihoods : hmm.get_genotyping_result()) {
				results->result.at(chromosome).at(index).combine(likelihoods);
				index += 1;
			}
		}
		// normalize the likelihoods after they have been combined
		for (size_t i = 0; i < results->result.at(chromosome).size(); ++i) {
			results->result.at(chromosome).at(i).normalize();
		}
	}
	// store runtime
	lock_guard<mutex> lock_result (results->result_mutex);
	if (results->runtimes.find(chromosome) == results->runtimes.end()) {
		results->runtimes.insert(pair<string,double>(chromosome, timer.get_total_time()));
	} else {
		results->runtimes[chromosome] += timer.get_total_time();
	}
}

bool ends_with (string const &full_string, string const ending) {
	if (full_string.size() >= ending.size()) {
		return (0 == full_string.compare(full_string.size() - ending.size(), ending.size(), ending));
	} else {
		return false;
	}
}

int main (int argc, char* argv[])
{
	Timer timer;
	double time_preprocessing;
	double time_kmer_counting;
	double time_path_sampling;
	double time_writing;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping and phasing based on kmer-counting and known haplotype sequences." << endl;
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
	long double effective_N = 0.00001L;
	long double regularization = 0.001L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	bool add_reference = true;
	size_t sampling_size = 10;
//	size_t sampling_times = 5;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
	argument_parser.add_optional_argument('n', "0.00001", "effective population size");
	argument_parser.add_flag_argument('g', "only run genotyping (Forward backward algorithm)");
	argument_parser.add_flag_argument('p', "only run phasing (Viterbi algorithm)");
	argument_parser.add_optional_argument('m', "0.001", "regularization constant for copynumber probabilities");
	argument_parser.add_flag_argument('c', "count all read kmers instead of only those located in graph.");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers.");
	argument_parser.add_flag_argument('d', "do not add reference as additional path.");
	argument_parser.add_optional_argument('a', "10", "sample subsets of paths of this size.");
//	argument_parser.add_optional_argument('b', "5", "how many times to sample paths.");

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
	effective_N = stold(argument_parser.get_argument('n'));
	regularization = stold(argument_parser.get_argument('m'));
	count_only_graph = !argument_parser.get_flag('c');
	ignore_imputed = argument_parser.get_flag('u');
	add_reference = !argument_parser.get_flag('d');
	sampling_size = stoi(argument_parser.get_argument('a'));
//	sampling_times = stoi(argument_parser.get_argument('b'));

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, kmersize, add_reference, sample_name);
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
//	cerr << "Compute kmer-coverage ..." << endl;
//	size_t kmer_coverage = read_kmer_counts.computeKmerCoverage(genome_kmers);
	size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
	cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

	// count kmers in allele + reference sequence
	cerr << "Count kmers in genome ..." << endl;
	KmerCounter* genomic_kmer_counts = new JellyfishCounter(segment_file, kmersize, nr_jellyfish_threads);

	// TODO: only for analysis
	struct rusage r_usage1;
	getrusage(RUSAGE_SELF, &r_usage1);
	cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

	// prepare output files
	if (! only_phasing) variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	if (! only_genotyping) variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	cerr << "Construct HMM and run core algorithm ..." << endl;
	time_kmer_counting = timer.get_interval_time();

	// prepare subsets of paths to run on
	size_t nr_paths = variant_reader.nr_of_paths();
	PathSampler path_sampler(nr_paths);
	vector<vector<size_t>> subsets;
	// TODO: determine how often to sample and how large each sample should be
//	path_sampler.select_multiple_subsets(subsets, sampling_size, sampling_times);
	path_sampler.partition_samples(subsets, sampling_size);

	for (auto s : subsets) {
		for (auto b : s) {
			cout << b << endl;
		}
		cout << "-----" << endl;
	}

	if (!only_phasing) cerr << "Sampled " << subsets.size() << " subset(s) of paths each of size " << sampling_size << " for genotyping." << endl;

	// for now, run phasing only once on largest set of paths that can still be handled.
	// in order to use all paths, an iterative stradegie should be considered
	vector<size_t> phasing_paths;
	size_t nr_phasing_paths = min((unsigned int) nr_paths, (unsigned int) 30);
	path_sampler.select_single_subset(phasing_paths, nr_phasing_paths);
	if (!only_genotyping) cerr << "Sampled " << phasing_paths.size() << " paths to be used for phasing." << endl;
	time_path_sampling = timer.get_interval_time();

	cerr << "Determine unique kmers ..." << endl;
	// determine number of cores to use
	size_t available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
	size_t nr_cores_uk = min(nr_core_threads, available_threads_uk);
	if (nr_cores_uk < nr_core_threads) {
		cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
	}

	// prepare UniqueKmers for each chromosome
	UniqueKmersMap unique_kmers_list;
	{
		// create thread pool with at most nr_chromosomes threads
		ThreadPool threadPool (nr_cores_uk);
		for (auto chromosome : chromosomes) {
			VariantReader* variants = &variant_reader;
			UniqueKmersMap* result = &unique_kmers_list;
			function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_kmer_counts, read_kmer_counts, variants, kmer_abundance_peak, regularization, result);
			threadPool.submit(f_unique_kmers);
		}
	}

	// TODO: only for analysis
	struct rusage r_usage2;
	getrusage(RUSAGE_SELF, &r_usage2);
	cerr << "#### Memory usage until now: " << (r_usage2.ru_maxrss / 1E6) << " GB ####" << endl;

	// destroy kmer counts as they are no longer needed
	delete genomic_kmer_counts;
	genomic_kmer_counts = nullptr;
	delete read_kmer_counts;
	read_kmer_counts = nullptr;

	// TODO: only for analysis
	struct rusage r_usage3;
	getrusage(RUSAGE_SELF, &r_usage3);
	cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;

	// determine max number of available threads for genotyping (at most one thread per chromosome and subsample possible)
	size_t available_threads = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size() * (unsigned int) subsets.size());
	if (nr_core_threads > available_threads) {
		cerr << "Warning: using " << available_threads << " for genotyping." << endl;
		nr_core_threads = available_threads;
	}

	// run genotyping
	Results results;
	{
		// create thread pool
		ThreadPool threadPool (nr_core_threads);
		for (auto chromosome : chromosomes) {
			vector<UniqueKmers*>* unique_kmers = &unique_kmers_list.unique_kmers[chromosome];
			Results* r = &results;
			// if requested, run phasing first
			if (!only_genotyping) {
				vector<size_t>* only_paths = &phasing_paths;
				function<void()> f_genotyping = bind(run_genotyping, chromosome, unique_kmers, false, true, effective_N, regularization, only_paths, r);
				threadPool.submit(f_genotyping);
			}

			if (!only_phasing) {
				// if requested, run genotying
				for (size_t s = 0; s < subsets.size(); ++s) {
					vector<size_t>* only_paths = &subsets[s];
					function<void()> f_genotyping = bind(run_genotyping, chromosome, unique_kmers, true, false, effective_N, regularization, only_paths, r);
					threadPool.submit(f_genotyping);
				}
			}
		}
	}

	timer.get_interval_time();

	// output VCF
	cerr << "Write results to VCF ..." << endl;
	if (!(only_genotyping && only_phasing)) assert (results.result.size() == chromosomes.size());
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
	cerr << "time spent counting kmers: \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent selecting paths: \t" << time_path_sampling << " sec" << endl;
	// output per chromosome time
	double time_hmm = time_writing;
	for (auto chromosome : chromosomes) {
		double time_chrom = results.runtimes[chromosome] + unique_kmers_list.runtimes[chromosome];
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << time_chrom << endl;
		time_hmm += time_chrom;
	}
	cerr << "total running time:\t" << time_preprocessing + time_kmer_counting + time_path_sampling + time_hmm + time_writing << " sec"<< endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
