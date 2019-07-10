#include <stdexcept>
#include <sstream>
#include "uniquekmers.hpp"

using namespace std;

UniqueKmers::UniqueKmers(size_t variant_id, size_t variant_position)
	:variant_id(variant_id),
	 variant_pos(variant_position),
	 current_index(0)
{}

size_t UniqueKmers::get_variant_index() {
	return this->variant_id;
}

size_t UniqueKmers::get_variant_position() {
	return this->variant_pos;
}

void UniqueKmers::insert_empty_path(size_t path) {
	this->paths[path] = KmerPath(); 
}

void UniqueKmers::insert_kmer(CopyNumber cn,  vector<size_t>& paths){
	size_t index = this->current_index;
	this->kmer_to_copynumber.push_back(cn);
	for (auto const& p: paths){
		this->paths[p].set_position(index);
	}
	current_index += 1;
}

bool UniqueKmers::kmer_on_path(size_t kmer_index, size_t path_index) const {
	if (kmer_index < this->current_index) {
		return (this->paths.at(path_index).get_position(kmer_index) > 0);
	} else {
		throw runtime_error("UniqueKmers::kmer_on_path: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

CopyNumberAssignment UniqueKmers::combine_paths(size_t path_id1, size_t path_id2) {
	return this->paths.at(path_id1) + this->paths.at(path_id2);
}

CopyNumber UniqueKmers::get_copynumber_of(size_t kmer_index) {
	if (kmer_index < this->current_index) {
		return this->kmer_to_copynumber[kmer_index];
	} else {
		throw runtime_error("UniqueKmers::get_copynumber_of: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

size_t UniqueKmers::size() const {
	return this->current_index;
}

void UniqueKmers::get_path_ids(std::vector<size_t>& result) {
	for (auto it = this->paths.begin(); it != this->paths.end(); ++it) {
		result.push_back(it->first);
	}
}

ostream& operator<< (ostream& stream, const UniqueKmers& uk) {
	stream << "UniqueKmers for variant: " << uk.variant_id << endl;
	for (size_t i = 0; i < uk.size(); ++i) {
		CopyNumber cn = uk.kmer_to_copynumber[i];
		stream << i << ": " << cn.get_probability_of(0) << " " << cn.get_probability_of(1) <<  " " << cn.get_probability_of(2) << endl;
	}
	stream << "paths:" << endl;
	for (auto it = uk.paths.begin(); it != uk.paths.end(); ++it) {
		stream << it->second.convert_to_string() << endl;
	}
	return stream;
}
