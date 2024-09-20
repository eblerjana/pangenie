#ifndef SAMPLED_PANEL
#define SAMPLED_PANEL

#include <vector>

class SampledPanel {
public:
	SampledPanel(std::vector<unsigned char> path_to_allele);
	unsigned char get_allele_on_path(size_t path) const;
	size_t get_nr_paths() const;
private:
	std::vector<unsigned char> path_to_allele;
};

#endif // SAMPLED_PANEL