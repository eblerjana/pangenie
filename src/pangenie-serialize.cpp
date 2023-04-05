#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <cereal/archives/binary.hpp>
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
	map<string, vector<shared_ptr<UniqueKmers>>> unique_kmers;
	map<string, double> runtimes;
};


struct Results {
	mutex result_mutex;
	map<string, shared_ptr<HMM>> result;
	map<string, double> runtimes;
};

void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, KmerCounter* read_kmer_counts, VariantReader* variant_reader, ProbabilityTable* probs, UniqueKmersMap* unique_kmers_map, size_t kmer_coverage) {
	Timer timer;
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, variant_reader, chromosome, kmer_coverage);
	std::vector<shared_ptr<UniqueKmers>> unique_kmers;
	kmer_computer.compute_unique_kmers(&unique_kmers, probs, true);
	// store the results
	{
		lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
		unique_kmers_map->unique_kmers.insert(pair<string, vector<shared_ptr<UniqueKmers>>> (chromosome, move(unique_kmers)));
	}
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes.insert(pair<string, double>(chromosome, timer.get_total_time()));
}

void run_genotyping(string chromosome, vector<shared_ptr<UniqueKmers>>* unique_kmers, ProbabilityTable* probs, bool only_genotyping, bool only_phasing, long double effective_N, vector<unsigned short>* only_paths, Results* results) {
	Timer timer;
	/* construct HMM and run genotyping/phasing. Genotyping is run without normalizing the final alpha*beta values.
	These values are first added up across different subsets of paths, and the resulting probabilities are normalized
	at the end. This is done so that genotyping runs on disjoint sets of paths are better comparable. */
	shared_ptr<HMM> hmm = shared_ptr<HMM>(new HMM(unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, false, effective_N, only_paths, false));
	// store the results
	{
		lock_guard<mutex> lock_result (results->result_mutex);
		// combine the new results to the already existing ones (if present)
		if (results->result.find(chromosome) == results->result.end()) {
			results->result.insert(pair<string, shared_ptr<HMM> > (chromosome, hmm) );
		} else {
			// combine newly computed likelihoods with already exisiting ones
			results->result.at(chromosome)->combine_likelihoods(*hmm);
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
	double time_preprocessing = 0.0;
	double time_kmer_counting = 0.0;
	double time_probabilities = 0.0;
	double time_unique_kmers = 0.0;
	double time_serialize = 0.0;
	double time_path_sampling = 0.0;
	double time_read_serialized = 0.0;
	double time_hmm = 0.0;
	double time_writing = 0.0;
	double time_total = 0.0;

	struct rusage rss_preprocessing;
	struct rusage rss_kmer_counting;
	struct rusage rss_probabilities;
	struct rusage rss_unique_kmers;
	struct rusage rss_serialize;
	struct rusage rss_path_sampling;
	struct rusage rss_read_serialized;
	struct rusage rss_hmm;
	struct rusage rss_total;


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
//	argument_parser.add_optional_argument('n', "0.00001", "effective population size");
	argument_parser.add_flag_argument('g', "run genotyping (Forward backward algorithm, default behaviour).");
	argument_parser.add_flag_argument('p', "run phasing (Viterbi algorithm). Experimental feature.");
//	argument_parser.add_optional_argument('m', "0.001", "regularization constant for copynumber probabilities");
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

//	effective_N = stold(argument_parser.get_argument('n'));
//	regularization = stold(argument_parser.get_argument('m'));
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


	// UniqueKmers for each chromosome
	UniqueKmersMap unique_kmers_list;
	ProbabilityTable probabilities;
	vector<string> chromosomes;
	unsigned short nr_paths = 0;

	{
	
		// read allele sequences and unitigs inbetween, write them into file
		cerr << "Determine allele sequences ..." << endl;
		VariantReader variant_reader (vcffile, reffile, kmersize, add_reference, sample_name);
				
		string segment_file = outname + "_path_segments.fasta";
		cerr << "Write path segments to file: " << segment_file << " ..." << endl;
		variant_reader.write_path_segments(segment_file);
		nr_paths = variant_reader.nr_of_paths();

		// determine chromosomes present in VCF
		variant_reader.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		getrusage(RUSAGE_SELF, &rss_preprocessing);
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
				read_kmer_counts = new JellyfishCounter(readfile, segment_file, kmersize, nr_jellyfish_threads, hash_size);
			} else {
				read_kmer_counts = new JellyfishCounter(readfile, kmersize, nr_jellyfish_threads, hash_size);
			}
		}

		size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
		cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

		// count kmers in allele + reference sequence
		cerr << "Count kmers in genome ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);

		getrusage(RUSAGE_SELF, &rss_kmer_counting);
		time_kmer_counting = timer.get_interval_time();

		cerr << "Serialize VariantReader object ..." << endl;
		{
  			ofstream os(outname + "_VariantReader.cereal", std::ios::binary);
  			cereal::BinaryOutputArchive archive( os );
			archive(variant_reader);
		}

		getrusage(RUSAGE_SELF, &rss_serialize);
		time_serialize = timer.get_interval_time();

		cerr << "Determine unique kmers ..." << endl;
		// determine number of cores to use
		size_t available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		size_t nr_cores_uk = min(nr_core_threads, available_threads_uk);
		if (nr_cores_uk < nr_core_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}

		// precompute probabilities
		probabilities = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);

		getrusage(RUSAGE_SELF, &rss_probabilities);
		time_probabilities = timer.get_interval_time();

		{
			// create thread pool with at most nr_chromosomes threads
			ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				VariantReader* variants = &variant_reader;
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				ProbabilityTable* probs = &probabilities;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, read_kmer_counts, variants, probs, result, kmer_abundance_peak);
				threadPool.submit(f_unique_kmers);
			}
		}

		delete read_kmer_counts;
		read_kmer_counts = nullptr;
	
		// compute total time spent determining unique kmers
		time_unique_kmers = 0.0;
		for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
			time_unique_kmers += it->second;
		}

		getrusage(RUSAGE_SELF, &rss_unique_kmers);
		timer.get_interval_time();
	}

	// prepare subsets of paths to run on
	// TODO: for too large panels, print waring
	if (nr_paths > 200) cerr << "Warning: panel is large and PanGenie might take a long time genotyping. Try reducing the panel size prior to genotyping." << endl;
	// handle case when sampling_size is not set
	if (sampling_size == 0) {
		if (nr_paths > 90) {
			sampling_size = 45;
		} else {
			sampling_size = nr_paths;		
		}
	}

	PathSampler path_sampler(nr_paths);
	vector<vector<unsigned short>> subsets;
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
	vector<unsigned short> phasing_paths;
	unsigned short nr_phasing_paths = min((unsigned short) nr_paths, (unsigned short) 30);
	path_sampler.select_single_subset(phasing_paths, nr_phasing_paths);
	if (!only_genotyping) cerr << "Sampled " << phasing_paths.size() << " paths to be used for phasing." << endl;

	getrusage(RUSAGE_SELF, &rss_path_sampling);
	time_path_sampling = timer.get_interval_time();

	cerr << "Construct HMM and run core algorithm ..." << endl;

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
			vector<shared_ptr<UniqueKmers>>* unique_kmers = &unique_kmers_list.unique_kmers[chromosome];
			ProbabilityTable* probs = &probabilities;
			Results* r = &results;
			// if requested, run phasing first
			if (!only_genotyping) {
				vector<unsigned short>* only_paths = &phasing_paths;
				function<void()> f_genotyping = bind(run_genotyping, chromosome, unique_kmers, probs, false, true, effective_N, only_paths, r);
				threadPool.submit(f_genotyping);
			}

			if (!only_phasing) {
				// if requested, run genotying
				for (size_t s = 0; s < subsets.size(); ++s){
					vector<unsigned short>* only_paths = &subsets[s];
					function<void()> f_genotyping = bind(run_genotyping, chromosome, unique_kmers, probs, true, false, effective_N, only_paths, r);
					threadPool.submit(f_genotyping);
				}
			}
		}
	}

	// in case genotyping was run, normalize the combined likelihoods
	if (!only_phasing){
		for (auto chromosome : chromosomes) {
			results.result.at(chromosome)->normalize();
		}
	}

	// compute total time spent genotyping
	for (auto it = results.runtimes.begin(); it != results.runtimes.end(); ++it) {
			time_hmm += it->second;
	}

	getrusage(RUSAGE_SELF, &rss_hmm);
	timer.get_interval_time();

	// read serialized VariantReader object
	VariantReader variant_reader;
	cerr << "Reading precomputed VariantReader ..." << endl; 
  	ifstream os(outname + "_VariantReader.cereal", std::ios::binary);
  	cereal::BinaryInputArchive archive( os );
	archive(variant_reader);

	getrusage(RUSAGE_SELF, &rss_read_serialized);
	time_read_serialized = timer.get_interval_time();


	// output VCF
	cerr << "Write results to VCF ..." << endl;

	// prepare output files
	if (! only_phasing) variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	if (! only_genotyping) variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	if (!(only_genotyping && only_phasing)) assert (results.result.size() == chromosomes.size());
	// write VCF
	for (auto it = results.result.begin(); it != results.result.end(); ++it) {
		if (!only_phasing) {
			// output genotyping results
			
			variant_reader.write_genotypes_of(it->first, it->second->get_ref_genotyping_result(), &unique_kmers_list.unique_kmers[it->first], ignore_imputed);
		}
		if (!only_genotyping) {
			// output phasing results
			variant_reader.write_phasing_of(it->first, it->second->get_ref_genotyping_result(), &unique_kmers_list.unique_kmers[it->first], ignore_imputed);
		}
	}

	if (! only_phasing) variant_reader.close_genotyping_outfile();
	if (! only_genotyping) variant_reader.close_phasing_outfile();

	getrusage(RUSAGE_SELF, &rss_total);
	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers (wallclock): \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent writing VariantReader to disk: \t" << time_serialize << " sec" << endl;
	cerr << "time spent pre-computing probabilities: \t" << time_probabilities << " sec" << endl;
	cerr << "time spent determining unique kmers: \t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent selecting paths: \t" << time_path_sampling << " sec" << endl;
	// output per chromosome time
	for (auto chromosome : chromosomes) {
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << results.runtimes[chromosome] << endl;
	}
	cerr << "time spent genotyping (total): \t" << time_hmm << endl;
	cerr << "time spent reading VariantReader from disk: \t" << time_read_serialized << " sec" << endl;
	cerr << "time spent writing output VCF: \t" << time_writing << " sec" << endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	cerr << endl;
	cerr << "Max RSS after reading input files: \t" << (rss_preprocessing.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers: \t" << (rss_kmer_counting.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after writing VariantReader to disk: \t" << (rss_serialize.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after pre-computing probabilities: \t" << (rss_probabilities.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after determining unique kmers: \t" << (rss_unique_kmers.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after selecting paths: \t" << (rss_path_sampling.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after genotyping: \t" << (rss_hmm.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after reading VariantReader from disk: \t" << (rss_read_serialized.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
