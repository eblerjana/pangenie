#ifndef COMMANDS_HPP
#define COMMANDS_HPP

#include <string>
#include <mutex>
#include <vector>
#include <memory>
#include <map>
#include "uniquekmers.hpp"

struct UniqueKmersMap {
	size_t kmersize;
	std::mutex kmers_mutex;
	std::map<std::string, std::vector<std::shared_ptr<UniqueKmers>>> unique_kmers;
	std::map<std::string, double> runtimes;
	bool add_reference;

	template <class Archive>
	void save(Archive& ar) const {
		ar(kmersize, unique_kmers, runtimes, add_reference);
	}

	template <class Archive>
	void load(Archive& ar) {
		ar(kmersize, unique_kmers, runtimes, add_reference);
	}
};

int run_single_command(std::string precomputed_prefix, std::string readfile, std::string reffile, std::string vcffile, size_t kmersize, std::string outname, std::string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, bool add_reference, size_t sampling_size, uint64_t hash_size, size_t panel_size, double recombrate, bool output_panel,  long double sampling_effective_N = 0.00001L, unsigned short allele_penalty = 10);

int run_index_command(std::string reffile, std::string vcffile, size_t kmersize, std::string outname, size_t nr_jellyfish_threads, bool add_reference, uint64_t hash_size);

int run_genotype_command(std::string precomputed_prefix, std::string readfile, std::string outname, std::string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, size_t sampling_size, uint64_t hash_size, size_t panel_size, double recombrate, bool output_panel, long double sampling_effective_N = 0.00001L, unsigned short allele_penalty = 10);

int run_sampling(std::string precomputed_prefix, std::string readfile, std::string outname, size_t nr_jellyfish_threads, size_t nr_core_threads, long double effective_N, long double regularization, bool count_only_graph, uint64_t hash_size, size_t panel_size, double recombrate, long double sampling_effective_N = 0.00001L, unsigned short allele_penalty = 10);


#endif // COMMANDS_HPP