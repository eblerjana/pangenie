#include <stdexcept>
#include <map>
#include "sampledpanel.hpp"


using namespace std;

SampledPanel::SampledPanel(vector<unsigned char> path_to_allele, size_t nr_unique_kmers)
	:path_to_allele(path_to_allele.size()),
	 unique_kmers(nr_unique_kmers)
{
	for (size_t i = 0; i < path_to_allele.size(); ++i) {
		this->path_to_allele[i] = (int) path_to_allele[i];
	}
}

SampledPanel::SampledPanel(vector<int> path_to_allele, size_t nr_unique_kmers)
	:path_to_allele(path_to_allele),
	 unique_kmers(nr_unique_kmers)
{}

int SampledPanel::get_allele_on_path(size_t path_id) const{
	if (path_id >= this->path_to_allele.size()) {
		throw runtime_error("SampledPanel:get_allele_on_path: index out of bounds.");
	}
	return this->path_to_allele.at(path_id);
}

size_t SampledPanel::get_nr_paths() const {
	return this->path_to_allele.size();
}

vector<int> SampledPanel::get_all_paths() const {
	return this->path_to_allele;
}

size_t SampledPanel::get_unique_kmers() const {
	return this->unique_kmers;
}

SampledPanel SampledPanel::get_specific_alleles(vector<unsigned char>& alleles) const {
	map<int, size_t> allele_to_idx;
	vector<int> updated_alleles(path_to_allele.size());
	for (size_t i = 0; i < alleles.size(); ++i) {
		allele_to_idx[alleles[i]] = i;
	}

	for (size_t i = 0; i < this->path_to_allele.size(); ++i) {
		int old_allele = path_to_allele[i];
		if (allele_to_idx.count(old_allele)) {
			updated_alleles[i] = allele_to_idx[old_allele];
		} else {
			updated_alleles[i] = -1;
		}
	}

	return SampledPanel(updated_alleles, this->unique_kmers);
}