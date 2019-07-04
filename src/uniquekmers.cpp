#include <stdexcept>
#include <sstream>
#include "uniquekmers.hpp"

using namespace std;

UniqueKmers::UniqueKmers(size_t variant_id, size_t variant_position)
	:variant_id(variant_id),
	 variant_pos(variant_position),
	 current_index(0)
{}

bool UniqueKmers::contains_kmer(string kmer) {
	map<string,size_t>::const_iterator it = this->kmer_to_index.find(kmer);
	if (it != this->kmer_to_index.end()){
		return true;
	} else {
		return false;
	}
}

bool UniqueKmers::kmer_on_path(string kmer, size_t path_id) {
	if (this->contains_kmer(kmer)){
		size_t index = this->kmer_to_index.at(kmer);
		return (this->paths.at(path_id).get_position(index) > 0 );
	} else {
		throw runtime_error("UniqueKmers::kmer_on_path: requested kmer: " + kmer + " does not exist.");
	}
}

CopyNumber UniqueKmers::get_copynumber_of(string kmer) {
	if (this->contains_kmer(kmer)){
		size_t index = this->kmer_to_index.at(kmer);
		return this->kmer_to_copynumber[index];
	} else {
		throw runtime_error("UniqueKmers::get_copynumber_of: requested kmer: " + kmer + " does not exist.");
	}
}

size_t UniqueKmers::get_variant_index() {
	return this->variant_id;
}

size_t UniqueKmers::get_variant_position() {
	return this->variant_pos;
}

void UniqueKmers::insert_empty_path(size_t path) {
	this->paths[path] = KmerPath(); 
}

void UniqueKmers::insert_kmer(string kmer, CopyNumber cn,  vector<size_t>& paths){
	size_t index = this->current_index;
	this->kmer_to_index.insert(pair<string,size_t>(kmer, index));
	this->kmer_to_copynumber.insert(pair<size_t,CopyNumber>(index,cn));
	for (auto const& p: paths){
		this->paths[p].set_position(index);
	}
	current_index += 1;
}

CopyNumberAssignment UniqueKmers::combine_paths(size_t path_id1, size_t path_id2) {
	return this->paths.at(path_id1) + this->paths.at(path_id2);
}

size_t UniqueKmers::size() {
	return this->current_index;
}

void UniqueKmers::get_path_ids(std::vector<size_t>& result) {
	for (auto it = this->paths.begin(); it != this->paths.end(); ++it) {
		result.push_back(it->first);
	}
}

ostream& operator<< (ostream& stream, const UniqueKmers& uk) {
	stream << "UniqueKmers for variant: " << uk.variant_id << endl; 
	for (auto it = uk.kmer_to_index.begin(); it != uk.kmer_to_index.end(); ++it) {
		CopyNumber cn = uk.kmer_to_copynumber.at(it->second);
		stream << it->second << ": " << it->first << " " << cn.get_probability_of(0) << " " << cn.get_probability_of(1) << " " << cn.get_probability_of(2) << endl;
	}
	return stream;
}
