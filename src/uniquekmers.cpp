#include <stdexcept>
#include <sstream>
#include "uniquekmers.hpp"
#include <math.h>

using namespace std;

UniqueKmers::UniqueKmers(size_t variant_id, size_t variant_position)
	:variant_id(variant_id),
	 variant_pos(variant_position),
	 current_index(0),
	 nr_paths(0),
	 max_path(0)
{}

size_t UniqueKmers::get_variant_index() {
	return this->variant_id;
}

size_t UniqueKmers::get_variant_position() {
	return this->variant_pos;
}

void UniqueKmers::insert_empty_allele(unsigned char allele_id) {
	this->alleles[allele_id] = KmerPath(); 
	this->allele_to_paths.insert(make_pair(allele_id,DynamicBitset()));
}

void UniqueKmers::insert_path(size_t path_id, unsigned char allele_id) {
	// unset bits that might be set already
	bool was_set = false;
	for (auto it = this->allele_to_paths.begin(); it != allele_to_paths.end(); ++it) {
		this->allele_to_paths[it->first].unset(path_id, was_set);
	}
	this->allele_to_paths[allele_id].set(path_id);
	if (!was_set) this->nr_paths += 1;
	this->max_path = max(this->max_path + 1, path_id + 1);
}

void UniqueKmers::insert_kmer(CopyNumber cn,  vector<unsigned char>& alleles){
	size_t index = this->current_index;
	this->kmer_to_copynumber.push_back(cn);
	for (auto const& a: alleles){
		this->alleles[a].set_position(index);
	}
	current_index += 1;
}

bool UniqueKmers::kmer_on_path(size_t kmer_index, size_t path_index) const {
	// check if path_id exists
	bool path_exists = false;
	unsigned char allele_id;
	for (auto it = this->allele_to_paths.begin(); it != this->allele_to_paths.end(); ++it) {
		if (it->second.is_set(path_index)) {
			path_exists = true;
			allele_id = it->first;
			break;
		}
	}
	if (!path_exists) {
		throw runtime_error("UniqueKmers::kmer_on_path: path_index " + to_string(path_index) + " does not exist.");
	}

	// check if kmer_index is valid and look up position
	if (kmer_index < this->current_index) {
		return (this->alleles.at(allele_id).get_position(kmer_index) > 0);
	} else {
		throw runtime_error("UniqueKmers::kmer_on_path: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

CopyNumberAssignment UniqueKmers::combine_paths(unsigned char allele_id1, unsigned char allele_id2) {
	return this->alleles.at(allele_id1) + this->alleles.at(allele_id2);
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

size_t UniqueKmers::get_nr_paths() const {
	return this->nr_paths;
}

void UniqueKmers::get_path_ids(vector<size_t>& p, vector<unsigned char>& a) {
	for (size_t path_id = 0; path_id < this->max_path; ++path_id) {
		for (auto it = this->allele_to_paths.begin(); it != this->allele_to_paths.end(); ++it) {
			// check which allele is covered (if any)
			if (it->second.is_set(path_id)) {
				p.push_back(path_id);
				a.push_back(it->first);
				break;
			}
		}
	}
}

void UniqueKmers::get_allele_ids(vector<unsigned char>& a) {
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		a.push_back(it->first);
	}
}

ostream& operator<< (ostream& stream, const UniqueKmers& uk) {
	stream << "UniqueKmers for variant: " << uk.variant_id << endl;
	for (size_t i = 0; i < uk.size(); ++i) {
		CopyNumber cn = uk.kmer_to_copynumber[i];
		stream << i << ": " << cn.get_probability_of(0) << " " << cn.get_probability_of(1) <<  " " << cn.get_probability_of(2) << endl;
	}
	stream << "alleles:" << endl;
	for (auto it = uk.alleles.begin(); it != uk.alleles.end(); ++it) {
		stream << it->first << "\t" << it->second.convert_to_string() << endl;
	}
	stream << "paths:" << endl;
	for (auto it = uk.allele_to_paths.begin(); it != uk.allele_to_paths.end(); ++it) {
		for (size_t path_id = 0; path_id < uk.max_path; ++path_id) {
			if (it->second.is_set(path_id)) stream << path_id << " covers allele " << it->first << endl;
		}
	}
	return stream;
}
