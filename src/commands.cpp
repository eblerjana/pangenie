#include "commands.hpp"

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
#include <zlib.h>
#include <cereal/archives/binary.hpp>
#include <cstdio>
#include "kmercounter.hpp"
#include "jellyfishreader.hpp"
#include "jellyfishcounter.hpp"
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "graph.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "pathsampler.hpp"
#include "kmerparser.hpp"
#include "graphbuilder.hpp"
#include "stepwiseuniquekmercomputer.hpp"
#include "uniquekmercomputer.hpp"
#include "threadpool.hpp"

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
	size_t kmersize;
	mutex kmers_mutex;
	map<string, vector<shared_ptr<UniqueKmers>>> unique_kmers;
	map<string, double> runtimes;

	template <class Archive>
	void save(Archive& ar) const {
		ar(kmersize, unique_kmers, runtimes);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(kmersize, unique_kmers, runtimes);
	}
};

struct Results {
	mutex result_mutex;
	map<string, vector<GenotypingResult>> result;
	map<string, double> runtimes;

	template <class Archive>
	void save(Archive& ar) const {
		ar(result, runtimes);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(result, runtimes);
	}
};


void fill_read_kmercounts_fasta(string chromosome, UniqueKmersMap* unique_kmers_map, shared_ptr<KmerCounter> read_kmer_counts, ProbabilityTable* probabilities, string outname, size_t kmer_coverage) {
	Timer timer;
	string filename = outname + "_" + chromosome + "_kmers.fa.gz";
	gzFile file = gzopen(filename.c_str(), "rb");
	if (!file) {
		throw runtime_error("fill_read_kmercounts_fasta: kmer file cannot be opened.");
	}

	const int buffer_size = 1024;
	char buffer[buffer_size];
	string line;
	size_t var_index = 0;
	size_t kmer_index = 0;
	bool is_flank = false;
	vector<string> kmers = {};
	vector<string> flanking_kmers = {};

    while (gzgets(file, buffer, buffer_size) != nullptr) {
        line += buffer;
        if (line.back() == '\n') {
			// remove newline character
			line.pop_back();

			if (line[0] == '>') {
				// FASTA header line. Assumes this format: ><u/f>_<kmer-index>_<variant-index>_<variant-position>
				vector<string> fields;
				parse(fields, line, '_');
				cout << line << endl;
				assert(fields.size() == 4);
				size_t read_var = atoi( fields[2].c_str() );
				kmer_index = atoi( fields[1].c_str() );
				if (read_var != var_index) {
					// new variant, store previously read kmer information
					size_t start = atoi(fields[3].c_str());
					assert(start == unique_kmers_map->unique_kmers[chromosome][read_var]->get_variant_position());

					{
						size_t kmers_used = 0;
						// add counts to UniqueKmers object
						for (size_t i = 0; i < kmers.size(); ++i) {
							assert (kmers_used <= 300);

							size_t count = read_kmer_counts->getKmerAbundance(kmers[i]);

							// determine probabilities
							CopyNumber cn = probabilities->get_probability(kmer_coverage, count);
							long double p_cn0 = cn.get_probability_of(0);
							long double p_cn1 = cn.get_probability_of(1);
							long double p_cn2 = cn.get_probability_of(2);

							// kmers with only 0 probabilities
							if (! ((p_cn0 > 0) || (p_cn1 > 0) || (p_cn2 > 0) )) cerr << "Warining: only zero probabilities for " << kmers[i] << " at " << chromosome << " " << start << endl; 

							kmers_used += 1;
							lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
							unique_kmers_map->unique_kmers[chromosome][var_index]->update_readcount(i, count);
							
						}
					}

					// determine local kmer coverage
					unsigned short local_coverage = compute_local_coverage(flanking_kmers, read_kmer_counts, kmer_coverage);

					lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
					unique_kmers_map->unique_kmers[chromosome][var_index]->set_coverage(local_coverage);

					// clear for next variant
					kmers = {};
					flanking_kmers = {};
					var_index = read_var;
					is_flank = false;
				} else {
					is_flank = (fields[0] == ">f");
				}
			} else {
				if (is_flank) {
					assert (kmer_index == flanking_kmers.size());
					flanking_kmers.push_back(line);
				} else {
					assert (kmer_index == kmers.size());
					kmers.push_back(line);
				}
			}
			line.clear();
        }
    }
	gzclose(file);
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes[chromosome] = timer.get_total_time();
}


void fill_read_kmercounts(string chromosome, UniqueKmersMap* unique_kmers_map, shared_ptr<KmerCounter> read_kmer_counts, ProbabilityTable* probabilities, string outname, size_t kmer_coverage) {
	Timer timer;
	string filename = outname + "_" + chromosome + "_kmers.tsv.gz";
	gzFile file = gzopen(filename.c_str(), "rb");
	if (!file) {
		throw runtime_error("fill_read_kmercounts: kmer file cannot be opened.");
	}

	const int buffer_size = 1024;
	char buffer[buffer_size];
	string line;
	size_t var_index = 0;
    while (gzgets(file, buffer, buffer_size) != nullptr) {
        line += buffer;
        if (line.back() == '\n') {

			// remove newline character
            line.pop_back();

			// read kmer information from file
			vector<string> kmers;
			vector<string> flanking_kmers;
			bool is_header = false;
			string chrom;
			size_t start;
			parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);

			// clear string for next line
            line.clear();

			if (is_header) continue; // header line
			assert(chrom == chromosome);
			assert(start == unique_kmers_map->unique_kmers[chromosome][var_index]->get_variant_position());

			{
				size_t kmers_used = 0;
				// add counts to UniqueKmers object
				for (size_t i = 0; i < kmers.size(); ++i) {
					assert (kmers_used <= 300);

					size_t count = read_kmer_counts->getKmerAbundance(kmers[i]);

					// determine probabilities
					CopyNumber cn = probabilities->get_probability(kmer_coverage, count);
					long double p_cn0 = cn.get_probability_of(0);
					long double p_cn1 = cn.get_probability_of(1);
					long double p_cn2 = cn.get_probability_of(2);

					if (!(p_cn0 > 0.0 || p_cn1 > 0.0 || p_cn2 > 0.0)) cerr << "Warining: only zero probabilities for " << kmers[i] << " at " << chrom << " " << start << endl; 

					kmers_used += 1;
					lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
					unique_kmers_map->unique_kmers[chromosome][var_index]->update_readcount(i, count);
				}
			}

			// determine local kmer coverage
			unsigned short local_coverage = compute_local_coverage(flanking_kmers, read_kmer_counts, kmer_coverage);

			lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
			unique_kmers_map->unique_kmers[chromosome][var_index]->set_coverage(local_coverage);
			var_index += 1;
        }
    }
	gzclose(file);
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes[chromosome] = timer.get_total_time();
}


void run_genotyping(string chromosome, vector<shared_ptr<UniqueKmers>>* unique_kmers, ProbabilityTable* probs, bool only_genotyping, bool only_phasing, long double effective_N, vector<unsigned short>* only_paths, Results* results) {
	Timer timer;
	/* construct HMM and run genotyping/phasing. Genotyping is run without normalizing the final alpha*beta values.
	These values are first added up across different subsets of paths, and the resulting probabilities are normalized
	at the end. This is done so that genotyping runs on disjoint sets of paths are better comparable. */
	HMM hmm(unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, false, effective_N, only_paths, false);

	// store the results
	{
		lock_guard<mutex> lock_result (results->result_mutex);
		// combine the new results to the already existing ones (if present)
		if (results->result.find(chromosome) == results->result.end()) {
			results->result.insert(pair<string, vector<GenotypingResult>> (chromosome, hmm.move_genotyping_result()));
		} else {
			// combine newly computed likelihoods with already exisiting ones
			size_t index = 0;
			vector<GenotypingResult> genotypes = hmm.move_genotyping_result();
			for (auto likelihoods : genotypes) {
				results->result.at(chromosome).at(index).combine(likelihoods);
				index += 1;
			}
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


void prepare_unique_kmers_stepwise(string chromosome, KmerCounter* genomic_kmer_counts, shared_ptr<Graph> graph, UniqueKmersMap* unique_kmers_map, string outname) {
	Timer timer;
	StepwiseUniqueKmerComputer kmer_computer(genomic_kmer_counts, graph);
	std::vector<shared_ptr<UniqueKmers>> unique_kmers;
	string filename = outname + "_" + chromosome + "_kmers.tsv.gz";
	kmer_computer.compute_unique_kmers(&unique_kmers, filename, true);
	// store the results
	{
		lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
		unique_kmers_map->unique_kmers.insert(pair<string, vector<shared_ptr<UniqueKmers>>> (chromosome, move(unique_kmers)));
	}
	// store runtime
	lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
	unique_kmers_map->runtimes.insert(pair<string, double>(chromosome, timer.get_total_time()));
}


// TODO: implement after adapting UniqueKmerComputer
void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, shared_ptr<KmerCounter> read_kmer_counts, shared_ptr<Graph> graph, ProbabilityTable* probs, UniqueKmersMap* unique_kmers_map, size_t kmer_coverage) {
	Timer timer;
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, graph, kmer_coverage);
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



int run_single_command(string precomputed_prefix, string readfile, string reffile, string vcffile, size_t kmersize, string outname, string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, bool add_reference, size_t sampling_size, uint64_t hash_size)
{

	Timer timer;
	double time_preprocessing = 0.0;
	double time_kmer_counting_graph = 0.0;
	double time_serialize_graph = 0.0;
	double time_unique_kmers = 0.0;
	double time_kmer_counting_reads = 0.0;
	double time_probabilities = 0.0;
	double time_path_sampling = 0.0;
	double time_hmm = 0.0;
	double time_writing = 0.0;
	double time_total = 0.0;

	struct rusage rss_preprocessing;
	struct rusage rss_kmer_counting_graph;
	struct rusage rss_serialize_graph;
	struct rusage rss_unique_kmers;
	struct rusage rss_kmer_counting_reads;
	struct rusage rss_probabilities;
	struct rusage rss_path_sampling;
	struct rusage rss_hmm;
	struct rusage rss_total;

	// check if input files exist and are uncompressed
	check_input_file(reffile);
	check_input_file(vcffile);
	check_input_file(readfile);

	UniqueKmersMap unique_kmers_list;
	ProbabilityTable probabilities;
	vector<string> chromosomes;
	Results results;
	string segment_file = outname + "_path_segments.fasta";
	size_t available_threads_uk;
	size_t nr_cores_uk;

	unique_kmers_list.kmersize = kmersize;

	/**
	*  1) Indexing step. Read variant information and determine unique kmers.
	*/

	{
		/** 
		*   Step 1: read variants, merge variants that are closer than kmersize apart,
		*   and write allele sequences and unitigs inbetween to a file.
		**/ 
		cerr << "Determine allele sequences ..." << endl;
		map<string, shared_ptr<Graph>> graph;
		GraphBuilder graph_builder (vcffile, reffile, graph, segment_file, kmersize, add_reference);

		// determine chromosomes present in VCF
		graph_builder.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		getrusage(RUSAGE_SELF, &rss_preprocessing);
		time_preprocessing = timer.get_interval_time();

		/**
		*  Step 2: count graph k-mers. Needed to determine unique k-mers in subsequent steps.
		**/ 
		cerr << "Count kmers in graph ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);


		getrusage(RUSAGE_SELF, &rss_kmer_counting_graph);
		time_kmer_counting_graph = timer.get_interval_time();


		/**
		* Step 3: count k-mers in reads
		**/

		shared_ptr<KmerCounter> read_kmer_counts = nullptr;
		// determine kmer copynumbers in reads
		if (readfile.substr(std::max(3, (int) readfile.size())-3) == std::string(".jf")) {
			cerr << "Read pre-computed read kmer counts ..." << endl;
			jellyfish::mer_dna::k(kmersize);
			read_kmer_counts = shared_ptr<JellyfishReader>(new JellyfishReader(readfile, kmersize));
		} else {
			cerr << "Count kmers in reads ..." << endl;

			if (count_only_graph) {
				read_kmer_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(readfile, {segment_file}, kmersize, nr_jellyfish_threads, hash_size));
			} else {
				read_kmer_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(readfile, kmersize, nr_jellyfish_threads, hash_size));
			}
		}


		/**
		* Step 4: Compute k-mer coverage and precompute probabilities.
		*/
		size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
		cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

		getrusage(RUSAGE_SELF, &rss_kmer_counting_reads);
		time_kmer_counting_reads = timer.get_interval_time();

		probabilities = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
		
		getrusage(RUSAGE_SELF, &rss_probabilities);
		time_probabilities = timer.get_interval_time();

		cerr << "Serialize Graph objects ..." << endl;
		for (auto chromosome : chromosomes) {
  			ofstream os(outname + "_" + chromosome + "_Graph.cereal", std::ios::binary);
  			cereal::BinaryOutputArchive archive( os );
			archive(*graph.at(chromosome));;
		}

		getrusage(RUSAGE_SELF, &rss_serialize_graph);
		time_serialize_graph = timer.get_interval_time();


		/**
		* Step 4: determine unique k-mers for each variant bubble 
		* This step will delete information from the GraphBuilder that will no longer be used.
		**/

		cerr << "Determine unique kmers ..." << endl;
		available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		nr_cores_uk = min(nr_jellyfish_threads, available_threads_uk);
		if (nr_cores_uk < nr_jellyfish_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}


		{
			// create thread pool with at most nr_chromosome threads
			ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				shared_ptr<Graph> graph_segment = graph.at(chromosome);
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				ProbabilityTable* probs = &probabilities;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, read_kmer_counts, graph_segment, probs, result, kmer_abundance_peak);
				threadPool.submit(f_unique_kmers);
			}
		}

		// determine the total runtime needed to compute unique kmers
		time_unique_kmers = 0.0;
		for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
			time_unique_kmers += it->second;
		}

		getrusage(RUSAGE_SELF, &rss_unique_kmers);
		timer.get_interval_time();	
	}


	/**
	*  2) Genotyping step. Construct a HMM and run the Forward-Backward algorithm to compute genotype likelihoods.
	*/

	{
		unsigned short nr_paths = 0;
		for (auto it = unique_kmers_list.unique_kmers.begin(); it != unique_kmers_list.unique_kmers.end(); ++it) {
			if (it->second.size() > 0) {
				nr_paths = it->second.at(0)->get_nr_paths();
			}
		}

		// TODO: for too large panels, print warning
		if (nr_paths > 500) cerr << "Warning: panel is large and PanGenie might take a long time genotyping. Try reducing the panel size prior to genotyping." << endl;
		// handle case when sampling_size is not set
		if (sampling_size == 0) {
			if (nr_paths > 220) {
				sampling_size = 110;
			} else {
				sampling_size = nr_paths;		
			}
		} else if (sampling_size > nr_paths) {
			// make sure that sampling size does not exceed nr_paths in panel
			sampling_size = nr_paths;
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
				for (size_t i = 0; i < results.result.at(chromosome).size(); ++i) {
					results.result.at(chromosome).at(i).normalize();
				}
			}
		}

		// compute total time spent genotyping
		for (auto it = results.runtimes.begin(); it != results.runtimes.end(); ++it) {
			time_hmm += it->second;
		}

		getrusage(RUSAGE_SELF, &rss_hmm);
		timer.get_interval_time();
	}

	// write the output VCF
	cerr << "Write results to VCF ..." << endl;
	if (!(only_genotyping && only_phasing)) assert (results.result.size() == chromosomes.size());
	bool write_header = true;
	for (auto it = results.result.begin(); it != results.result.end(); ++it) {
		// read serialized Graph object corresponding to current chromosome
		Graph graph;
        string graph_filename = precomputed_prefix + "_" + it->first + "_Graph.cereal";
		cerr << "Reading precomputed Graph for chromosome " << it->first << " ..." <<  " from " << graph_filename << endl;
		ifstream os(graph_filename, std::ios::binary);
		cereal::BinaryInputArchive archive( os );
		archive(graph);

		cerr << "Writing results for chromosome " << it->first << " ..." << endl;
		if (!only_phasing) {
			// output genotyping results
			graph.write_genotypes(outname + "_genotyping.vcf", it->second, write_header, sample_name, ignore_imputed);
		}
		if (!only_genotyping) {
			// output phasing results
			graph.write_phasing(outname + "_phasing.vcf", it->second, write_header, sample_name, ignore_imputed);
		}
		// write header only for first chromosome
		write_header = false;

		// remove file
		remove(graph_filename.c_str());
	}

	getrusage(RUSAGE_SELF, &rss_total);
	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();


	cerr << endl << "############### Summary ###############" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers in genome (wallclock): \t" << time_kmer_counting_graph << " sec" << endl;
	cerr << "time spent counting kmers in reads (wallclock): \t" << time_kmer_counting_reads << " sec" << endl;
	cerr << "time spent pre-computing probabilities: \t" << time_probabilities << " sec" << endl;
	cerr << "time spent writing Graph objects to disk: \t" << time_serialize_graph << " sec" << endl;
	cerr << "time spent determining unique kmers: \t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent selecting paths: \t" << time_path_sampling << " sec" << endl;
	// output per chromosome time
	for (auto chromosome : chromosomes) {
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << results.runtimes[chromosome] << endl;
	}
	cerr << "time spent genotyping (total): \t" << time_hmm << " sec" << endl;
	cerr << "time spent writing output VCF: \t" << time_writing << " sec" << endl;
	cerr << "total wallclock time PanGenie: " << time_total  << " sec" << endl;

	cerr << endl;
	cerr << "Max RSS after reading input files: \t" << (rss_preprocessing.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers in genome: \t" << (rss_kmer_counting_graph.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers in reads: \t" << (rss_kmer_counting_reads.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after pre-computing probabilities: \t" << (rss_probabilities.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after determining unique kmers: \t" << (rss_unique_kmers.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after selecting paths: \t" << (rss_path_sampling.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after genotyping: \t" << (rss_hmm.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
    cerr << "#######################################" << endl << endl;

	return 0;
}


int run_index_command(string reffile, string vcffile, size_t kmersize, string outname, size_t nr_jellyfish_threads, bool add_reference, uint64_t hash_size)
{

	Timer timer;
	double time_preprocessing = 0.0;
	double time_kmer_counting = 0.0;
	double time_serialize_graph = 0.0;
	double time_unique_kmers = 0.0;
	double time_serialize = 0.0;
	double time_total = 0.0;

	struct rusage rss_preprocessing;
	struct rusage rss_kmer_counting;
	struct rusage rss_serialize_graph;
	struct rusage rss_unique_kmers;
	struct rusage rss_total;

	// check if input files exist and are uncompressed
	check_input_file(reffile);
	check_input_file(vcffile);

	UniqueKmersMap unique_kmers_list;
	vector<string> chromosomes;
	string segment_file = outname + "_path_segments.fasta";
	size_t available_threads_uk;
	size_t nr_cores_uk;

	unique_kmers_list.kmersize = kmersize;

	/**
	*  1) Indexing step. Read variant information and determine unique kmers.
	*/

	{
		/** 
		*   Step 1: read variants, merge variants that are closer than kmersize apart,
		*   and write allele sequences and unitigs inbetween to a file.
		**/ 
		cerr << "Determine allele sequences ..." << endl;
		map<string, shared_ptr<Graph>> graph;
		GraphBuilder graph_builder (vcffile, reffile, graph, segment_file, kmersize, add_reference);

		// determine chromosomes present in VCF
		graph_builder.get_chromosomes(&chromosomes);
		cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

		getrusage(RUSAGE_SELF, &rss_preprocessing);
		time_preprocessing = timer.get_interval_time();

		/**
		*  Step 2: count graph k-mers. Needed to determine unique k-mers in subsequent steps.
		**/ 
		cerr << "Count kmers in graph ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);


		getrusage(RUSAGE_SELF, &rss_kmer_counting);
		time_kmer_counting = timer.get_interval_time();

		cerr << "Serialize Graph objects ..." << endl;
		for (auto chromosome : chromosomes) {
  			ofstream os(outname + "_" + chromosome + "_Graph.cereal", std::ios::binary);
  			cereal::BinaryOutputArchive archive( os );
			archive(*graph.at(chromosome));;
		}

		getrusage(RUSAGE_SELF, &rss_serialize_graph);
		time_serialize_graph = timer.get_interval_time();


		/**
		* Step 3: determine unique k-mers for each variant bubble and prepare datastructure storing
		* that information. It will later be used for the genotyping step. 
		* This step will delete information from the GraphBuilder that will no longer be used.
		**/

		cerr << "Determine unique kmers ..." << endl;
		available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		nr_cores_uk = min(nr_jellyfish_threads, available_threads_uk);
		if (nr_cores_uk < nr_jellyfish_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}


		{
			// create thread pool with at most nr_chromosome threads
			ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				shared_ptr<Graph> graph_segment = graph.at(chromosome);
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers_stepwise, chromosome, genomic_counts, graph_segment, result, outname);
				threadPool.submit(f_unique_kmers);
			}
		}

		// determine the total runtime needed to compute unique kmers
		time_unique_kmers = 0.0;
		for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
			time_unique_kmers += it->second;
		}

		getrusage(RUSAGE_SELF, &rss_unique_kmers);
		timer.get_interval_time();	
	}

	// serialization of UniqueKmersMap object
	cerr << "Storing unique kmer information ..." << endl;
	{
  		ofstream os(outname + "_UniqueKmersMap.cereal", std::ios::binary);
  		cereal::BinaryOutputArchive archive( os );
		archive(unique_kmers_list);
	}

	getrusage(RUSAGE_SELF, &rss_total);
	time_serialize = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary PanGenie-index ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers in genome (wallclock): \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent writing Graph objects to disk: \t" << time_serialize_graph << " sec" << endl;
	cerr << "time spent determining unique kmers: \t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent writing UniqueKmersMap to disk: \t" << time_serialize << " sec" << endl;
	cerr << "total wallclock time PanGenie-index: " << time_total  << " sec" << endl;

	cerr << endl;
	cerr << "Max RSS after reading input files: \t" << (rss_preprocessing.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers in genome: \t" << (rss_kmer_counting.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after determining unique kmers: \t" << (rss_unique_kmers.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
    cerr << "####################################" << endl << endl;
	return 0;

}

int run_genotype_command(string precomputed_prefix, string readfile, string outname, string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, size_t sampling_size, uint64_t hash_size)
{

	Timer timer;
	double time_read_serialized = 0.0;
	double time_unique_kmers = 0.0;
	double time_kmer_counting = 0.0;
	double time_probabilities = 0.0;
	double time_path_sampling = 0.0;
	double time_hmm = 0.0;
	double time_writing = 0.0;
	double time_total = 0.0;

	struct rusage rss_read_serialized;
	struct rusage rss_unique_kmers;
	struct rusage rss_kmer_counting;
	struct rusage rss_probabilities;
	struct rusage rss_path_sampling;
	struct rusage rss_hmm;
	struct rusage rss_total;

	// check if input files exist and are uncompressed
	check_input_file(readfile);

	vector<string> chromosomes;
	Results results;

	{
		UniqueKmersMap unique_kmers_list;
		ProbabilityTable probabilities;
		string segment_file = precomputed_prefix + "_path_segments.fasta";
        check_input_file(segment_file);
		size_t available_threads_uk;
		size_t nr_cores_uk;
		unsigned short nr_paths = 0;


		// re-construct UniqueKmersMap + chromosomes from input file (-f)
		string unique_kmers_archive = precomputed_prefix + "_UniqueKmersMap.cereal";
        check_input_file(unique_kmers_archive);
		cerr << "Reading precomputed UniqueKmersMap from " << unique_kmers_archive << " ..." << endl;
		ifstream is(unique_kmers_archive, std::ios::binary);
		cereal::BinaryInputArchive archive_is( is );
		archive_is(unique_kmers_list);

		// check if there are any variants
		size_t variants_read = 0;
		for (auto it = unique_kmers_list.unique_kmers.begin(); it != unique_kmers_list.unique_kmers.end(); ++it) {
			chromosomes.push_back(it->first);
			if (it->second.size() > 0) {
				nr_paths = it->second.at(0)->get_nr_paths();
				variants_read += it->second.size();
			}
		}

		cerr << "Read " << variants_read << " variants from provided UniqueKmersMap archive." << endl;

		// if no variants present, nothing to be done. Exit program.
		if (variants_read == 0) return 0;

		// there must be paths given.
		if (nr_paths == 0) {
			throw runtime_error("PanGenie-index: no haplotype paths given.");
		}

		getrusage(RUSAGE_SELF, &rss_read_serialized);
		time_read_serialized = timer.get_interval_time();

		/**
		*  2) K-mer counting in sequencing reads
		*/

		{
			/**
			* Step 1: Count kmers in the sequencing reads of the sample using Jellyfish
			* or read already computed counts from .jf file.
			*/

			size_t kmersize = unique_kmers_list.kmersize;

			shared_ptr<KmerCounter> read_kmer_counts = nullptr;
			// determine kmer copynumbers in reads
			if (readfile.substr(std::max(3, (int) readfile.size())-3) == std::string(".jf")) {
				cerr << "Read pre-computed read kmer counts ..." << endl;
				jellyfish::mer_dna::k(kmersize);
				read_kmer_counts = shared_ptr<JellyfishReader>(new JellyfishReader(readfile, kmersize));
			} else {
				cerr << "Count kmers in reads ..." << endl;

				if (count_only_graph) {
					read_kmer_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(readfile, {segment_file}, kmersize, nr_jellyfish_threads, hash_size));
				} else {
					read_kmer_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(readfile, kmersize, nr_jellyfish_threads, hash_size));
				}
			}


			/**
			* Step 2: Compute k-mer coverage and precompute probabilities.
			*/
			size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
			cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

			getrusage(RUSAGE_SELF, &rss_kmer_counting);
			time_kmer_counting = timer.get_interval_time();

			probabilities = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
		
			getrusage(RUSAGE_SELF, &rss_probabilities);
			time_probabilities = timer.get_interval_time();


			/**
			* Step 3: Fill the UniqueKmers object with sample-specific kmer counts.
			* Also, estimate the local kmer coverage from nearby kmers.
			*/

			cerr << "Determine read k-mer counts for unique kmers ..." << endl;
			available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
			nr_cores_uk = min(nr_core_threads, available_threads_uk);
			if (nr_cores_uk < nr_core_threads) {
				cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
			}
			
			{
				ThreadPool threadPool (nr_cores_uk);
				for (auto chromosome : chromosomes) {
					UniqueKmersMap* unique_kmers = &unique_kmers_list;
					ProbabilityTable* probs = &probabilities;
					function<void()> f_fill_readkmers = bind(fill_read_kmercounts, chromosome, unique_kmers, read_kmer_counts, probs, precomputed_prefix, kmer_abundance_peak);
					threadPool.submit(f_fill_readkmers);
				}
			}

			// determine the total runtime needed to update kmer information
			for (auto it = unique_kmers_list.runtimes.begin(); it != unique_kmers_list.runtimes.end(); ++it) {
				time_unique_kmers += it->second;
			}

			getrusage(RUSAGE_SELF, &rss_unique_kmers);
			timer.get_interval_time();

		}


		/**
		* 3) Genotyping. Construct a HMM and run the Forward-Backward algorithm to compute genotype likelihoods.
		*/

		{
			// TODO: for too large panels, print warning
			if (nr_paths > 500) cerr << "Warning: panel is large and PanGenie might take a long time genotyping. Try reducing the panel size prior to genotyping." << endl;
			// handle case when sampling_size is not set
			if (sampling_size == 0) {
				if (nr_paths > 220) {
					sampling_size = 110;
				} else {
					sampling_size = nr_paths;		
				}
			} else if (sampling_size > nr_paths) {
				// make sure that sampling size does not exceed nr_paths in panel
				sampling_size = nr_paths;
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
					for (size_t i = 0; i < results.result.at(chromosome).size(); ++i) {
						results.result.at(chromosome).at(i).normalize();
				}
				}
			}

			// compute total time spent genotyping
			for (auto it = results.runtimes.begin(); it != results.runtimes.end(); ++it) {
				time_hmm += it->second;
			}

			getrusage(RUSAGE_SELF, &rss_hmm);
			timer.get_interval_time();
		}
	}

	// write the output VCF
	cerr << "Write results to VCF ..." << endl;
	if (!(only_genotyping && only_phasing)) assert (results.result.size() == chromosomes.size());
	bool write_header = true;
	for (auto it = results.result.begin(); it != results.result.end(); ++it) {
		// read serialized Graph object corresponding to current chromosome
		Graph graph;
        string graph_filename = precomputed_prefix + "_" + it->first + "_Graph.cereal";
		cerr << "Reading precomputed Graph for chromosome " << it->first << " ..." <<  " from " << graph_filename << endl;
		ifstream os(graph_filename, std::ios::binary);
		cereal::BinaryInputArchive archive( os );
		archive(graph);

		cerr << "Writing results for chromosome " << it->first << " ..." << endl;
		if (!only_phasing) {
			// output genotyping results
			graph.write_genotypes(outname + "_genotyping.vcf", it->second, write_header, sample_name, ignore_imputed);
		}
		if (!only_genotyping) {
			// output phasing results
			graph.write_phasing(outname + "_phasing.vcf", it->second, write_header, sample_name, ignore_imputed);
		}
		// write header only for first chromosome
		write_header = false;
	}

	getrusage(RUSAGE_SELF, &rss_total);
	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary PanGenie-genotype ######" << endl;
	// output times
	cerr << "time spent reading UniqueKmersMap from disk: \t" << time_read_serialized << " sec" << endl;
	cerr << "time spent counting kmers in reads (wallclock): \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent pre-computing probabilities: \t" << time_probabilities << " sec" << endl;
	cerr << "time spent updating unique kmers: \t" << time_unique_kmers << " sec" << endl;
	cerr << "time spent selecting paths: \t" << time_path_sampling << " sec" << endl;
	// output per chromosome time
	for (auto chromosome : chromosomes) {
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << results.runtimes[chromosome] << endl;
	}
	cerr << "time spent genotyping (total): \t" << time_hmm << " sec" << endl;

	cerr << "time spent writing output VCF: \t" << time_writing << " sec" << endl;
	cerr << "total wallclock time PanGenie-genotype: " << time_total  << " sec" << endl;

	cerr << endl;
	cerr << "Max RSS after reading UniqueKmersMap from disk: \t" << (rss_read_serialized.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after counting kmers in reads: \t" << (rss_kmer_counting.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after pre-computing probabilities: \t" << (rss_probabilities.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after updating unique kmers: \t" << (rss_unique_kmers.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after selecting paths: \t" << (rss_path_sampling.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS after genotyping: \t" << (rss_hmm.ru_maxrss / 1E6) << " GB" << endl;
	cerr << "Max RSS: \t" << (rss_total.ru_maxrss / 1E6) << " GB" << endl;
    cerr << "#######################################" << endl << endl;
	return 0;

}
