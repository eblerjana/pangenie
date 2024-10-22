#include <stdexcept>
#include <sstream>
#include "uniquekmers.hpp"

using namespace std;

UniqueKmers::UniqueKmers(size_t variant_position, vector<unsigned char>& alleles)
	:variant_pos(variant_position),
	 current_index(0),
	 path_to_allele(alleles.size()),
	 local_coverage(0)
{
	for (size_t i = 0; i < alleles.size(); ++i) {
		unsigned char a = alleles[i];
		this->path_to_allele[i] = a;
		this->alleles[a] = AlleleInfo();
	}
}

size_t UniqueKmers::get_variant_position() {
	return this->variant_pos;
}


// NOTE: hardcoded threshold of 2 used here!
void UniqueKmers::insert_kmer(unsigned short readcount,  vector<unsigned char>& alleles){
	size_t index = this->current_index;
	this->kmer_to_count.push_back(readcount);
	for (auto const& a: alleles){
		this->alleles[a].kmer_path.set_position(index);
	}
	current_index += 1;
}

bool UniqueKmers::kmer_on_path(size_t kmer_index, size_t path_index) const {
	// check if path_id exists
	if (path_index >= this->path_to_allele.size()) {
		throw runtime_error("UniqueKmers::kmer_on_path: path_index " + to_string(path_index) + " does not exist.");
	}

	// check if kmer_index is valid and look up position
	if (kmer_index < this->current_index) {
		unsigned char allele_id = this->path_to_allele.at(path_index);
		return (this->alleles.at(allele_id).kmer_path.get_position(kmer_index) > 0);
	} else {
		throw runtime_error("UniqueKmers::kmer_on_path: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

unsigned short UniqueKmers::get_readcount_of(size_t kmer_index) {
	if (kmer_index < this->current_index) {
		return this->kmer_to_count[kmer_index];
	} else {
		throw runtime_error("UniqueKmers::get_readcount_of: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

void UniqueKmers::update_readcount(size_t kmer_index, unsigned short new_count) {
	if (kmer_index <  this->current_index) {
		this->kmer_to_count[kmer_index] = new_count;
	} else {
		throw runtime_error("UniqueKmers::update_readcount: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

size_t UniqueKmers::size() const {
	return this->current_index;
}

unsigned short UniqueKmers::get_nr_paths() const {
	return this->path_to_allele.size();
}

void UniqueKmers::get_path_ids(vector<unsigned short>& p, vector<unsigned char>& a, vector<unsigned short>* only_include) {
	if (only_include != nullptr) {
		// only return paths that are also contained in only_include
		for (auto p_it = only_include->begin(); p_it != only_include->end(); ++p_it) {
			// check if path is in path_to_allele
			if (*p_it < this->path_to_allele.size()) {
				p.push_back(*p_it);
				a.push_back(this->path_to_allele.at(*p_it));
			}
		}
	} else {
		// return all paths and corresponding alleles
		for (size_t i = 0; i < this->path_to_allele.size(); ++i) {
			p.push_back(i);
			a.push_back(this->path_to_allele.at(i));
		}
	}
}

void UniqueKmers::get_allele_ids(vector<unsigned char>& a) {
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		a.push_back(it->first);
	}
}

void UniqueKmers::get_defined_allele_ids(std::vector<unsigned char>& a) {
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		if (!it->second.is_undefined) a.push_back(it->first);
	}
}

ostream& operator<< (ostream& stream, const UniqueKmers& uk) {
	stream << "UniqueKmers for variant: " << uk.variant_pos << endl;
	for (size_t i = 0; i < uk.size(); ++i) {
		stream << i << ": " << uk.kmer_to_count[i] << endl;
	}
	stream << "alleles:" << endl;
	for (auto it = uk.alleles.begin(); it != uk.alleles.end(); ++it) {
		stream << (unsigned int) it->first << "\t" << it->second.kmer_path.convert_to_string() << endl;
	}

	stream << "undefined alleles:" << endl;
	for (auto it = uk.alleles.begin(); it != uk.alleles.end(); ++it) {
		if (uk.is_undefined_allele(it->first)) stream << (unsigned int) it->first << endl;
	}

	stream << "paths:" << endl;
	for (size_t  i = 0; i < uk.path_to_allele.size(); ++i) {
		stream << i << " covers allele " << (unsigned int) uk.path_to_allele.at(i) << endl;
	}

	stream << "local coverage: " << uk.local_coverage << endl;
	return stream;
}

void UniqueKmers::set_coverage(unsigned short local_coverage) {
	this->local_coverage = local_coverage;
}

unsigned short UniqueKmers::get_coverage() const {
	return this->local_coverage;
}

map<unsigned char, int> UniqueKmers::kmers_on_alleles () const {
	map<unsigned char, int> result;
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		result[it->first] = alleles.at(it->first).kmer_path.nr_kmers();
	}
	return result;
}


unsigned short UniqueKmers::kmers_on_allele(unsigned char allele_id) const {
	return alleles.at(allele_id).kmer_path.nr_kmers();
}


unsigned short UniqueKmers::present_kmers_on_allele(unsigned char allele_id) const {
	unsigned short result = 0;
	for (size_t i = 0; i < this->kmer_to_count.size(); ++i) {
		if (kmer_to_count[i] < 3) continue;
		if (alleles.at(allele_id).kmer_path.get_position(i) > 0) result += 1;
	}
	return result;
}

float UniqueKmers::fraction_present_kmers_on_allele(unsigned char allele_id) const {
	unsigned short total = this->kmers_on_allele(allele_id);
	if (total > 0) return this->present_kmers_on_allele(allele_id) / (float) total;
	return 1.0;
}

bool UniqueKmers::is_undefined_allele (unsigned char allele_id) const {
	// check if allele id exists
	auto it = this->alleles.find(allele_id);
	if (it != this->alleles.end()) {
		return it->second.is_undefined;
	} else {
		return false;
	}
}

void UniqueKmers::set_undefined_allele (unsigned char allele_id) {
	auto it = this->alleles.find(allele_id);
	if (it == this->alleles.end()) {
		throw runtime_error("UniqueKmers::set_undefined_allele: allele_id " + to_string(allele_id) + " does not exist.");
	}
	this->alleles[allele_id].is_undefined = true;
}

unsigned char UniqueKmers::get_allele(unsigned short path_id) const {
	if (path_id >= this->path_to_allele.size()) {
		throw runtime_error("UniqueKmers:get_allele: index out of bounds.");
	}
	return this->path_to_allele[path_id];
}

void UniqueKmers::update_paths(vector<unsigned short>& path_ids) {
	size_t nr_paths = path_ids.size();
	vector<unsigned char> updated_path_to_allele(nr_paths);
	map<unsigned char, AlleleInfo> updated_alleles;
	vector<unsigned char> undefined_alleles;
	for (size_t i = 0; i < path_ids.size(); ++i) {
		unsigned char allele = this->get_allele(path_ids[i]);
		updated_path_to_allele[i] = allele;
		updated_alleles[allele] = this->alleles[allele];
	}

	// update the KmerPath objects and the kmer counts
	map<size_t, vector<unsigned char>> kmer_to_alleles;
	for (auto it = updated_alleles.begin(); it != updated_alleles.end(); ++it) {
		for (size_t k = 0; k < this->size(); ++k) {
			if (it->second.kmer_path.get_position(k)) {
				kmer_to_alleles[k].push_back(it->first);
			}
		}
		if (it->second.is_undefined) undefined_alleles.push_back(it->first);
	}
	this->path_to_allele = updated_path_to_allele;
	this->alleles.clear();
	for (auto a : updated_path_to_allele) {
		this->alleles[a] = AlleleInfo();
	}

	vector<unsigned short> old_counts = this->kmer_to_count;
	this->kmer_to_count.clear();
	this->current_index = 0;
	// set undefined alleles
	for (auto a : undefined_alleles) this->set_undefined_allele(a);

	for (auto it = kmer_to_alleles.begin(); it != kmer_to_alleles.end(); ++it) {
		this->insert_kmer(old_counts[it->first], it->second);
	}
}