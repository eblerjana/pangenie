#ifndef COMMANDS_HPP
#define COMMANDS_HPP

#include <string>

int run_single_command(std::string precomputed_prefix, std::string readfile, std::string reffile, std::string vcffile, size_t kmersize, std::string outname, std::string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, bool add_reference, size_t sampling_size, uint64_t hash_size);

int run_index_command(std::string reffile, std::string vcffile, size_t kmersize, std::string outname, size_t nr_jellyfish_threads, bool add_reference, uint64_t hash_size);

int run_genotype_command(std::string precomputed_prefix, std::string readfile, std::string outname, std::string sample_name, size_t nr_jellyfish_threads, size_t nr_core_threads, bool only_genotyping, bool only_phasing, long double effective_N, long double regularization, bool count_only_graph, bool ignore_imputed, size_t sampling_size, uint64_t hash_size);

#endif // COMMANDS_HPP