#include <stdexcept>
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