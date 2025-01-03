#ifndef SAMPLED_PANEL
#define SAMPLED_PANEL

#include <vector>

class SampledPanel {
public:
	SampledPanel(std::vector<unsigned char> path_to_allele, size_t nr_unique_kmers);
	SampledPanel(std::vector<int> path_to_allele, size_t nr_unique_kmers);
	int get_allele_on_path(size_t path) const;
	std::vector<int> get_all_paths() const;
	size_t get_unique_kmers() const;
	size_t get_nr_paths() const;
	SampledPanel get_specific_alleles(std::vector<unsigned char>& alleles) const;
private:
	std::vector<int> path_to_allele;
	size_t unique_kmers;
};

#endif // SAMPLED_PANEL