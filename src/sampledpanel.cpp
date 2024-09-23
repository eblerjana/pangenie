#include <stdexcept>
#include <map>
#include "sampledpanel.hpp"


using namespace std;

SampledPanel::SampledPanel(vector<unsigned char> path_to_allele)
	:path_to_allele(path_to_allele)
{}

unsigned char SampledPanel::get_allele_on_path(size_t path_id) const{
	if (path_id >= this->path_to_allele.size()) {
		throw runtime_error("SampledPanel:get_allele_on_path: index out of bounds.");
	}
	return this->path_to_allele.at(path_id);
}

size_t SampledPanel::get_nr_paths() const {
	return this->path_to_allele.size();
}

vector<unsigned char> SampledPanel::get_all_paths() const {
	return this->path_to_allele;
}

SampledPanel SampledPanel::get_specific_alleles(vector<unsigned char>& alleles) const {
	map<unsigned char, size_t> allele_to_idx;
	vector<unsigned char> updated_alleles(alleles.size());
	for (size_t i = 0; i < alleles.size(); ++i) {
		allele_to_idx[alleles[i]] = i;
	}

	for (size_t i = 0; i < this->path_to_allele.size(); ++i) {
		unsigned char old_allele = path_to_allele[i];
		if (allele_to_idx.count(old_allele)) updated_alleles[allele_to_idx[old_allele]];
	}

	return SampledPanel(updated_alleles);
}