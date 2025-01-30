#include <stdexcept>
#include <sstream>
#include <cassert>
#include "biallelicuniquekmers.hpp"

using namespace std;

BiallelicUniqueKmers::BiallelicUniqueKmers(size_t variant_position, vector<unsigned short>& alleles)
	:UniqueKmers(variant_position, 0),
	 current_index(0),
	 path_to_allele(alleles.size())
{

	for (size_t i = 0; i < alleles.size(); ++i) {
		unsigned short a = alleles[i];
		if ((a != 0) && (a != 1)) {
			throw runtime_error("BiallelicUniqueKmers::BiallelicUniqueKmers: provided alleles need to be either 0 or 1 (biallelic).");
		}
		bool allele = (bool) a;
		this->path_to_allele[i] = allele;
		this->alleles[allele] = AlleleInfo();
	}
}


void BiallelicUniqueKmers::insert_kmer(unsigned short readcount,  vector<unsigned short>& alleles){
	size_t index = this->current_index;
	this->kmer_to_count.push_back(readcount);
	for (auto const& a: alleles){
		if ((a != 1) && (a != 0)) {
			cout << "A IS THIS: " << a << endl;
			throw runtime_error("BiallelicUniqueKmers::insert_kmer: provided alleles need to be either 0 or 1 (biallelic)");
		}
		this->alleles[a].kmer_path.set_position(index);
	}
	current_index += 1;
}

bool BiallelicUniqueKmers::kmer_on_path(size_t kmer_index, size_t path_index) const {
	// check if path_id exists
	if (path_index >= this->path_to_allele.size()) {
		throw runtime_error("BiallelicUniqueKmers::kmer_on_path: path_index " + to_string(path_index) + " does not exist.");
	}

	// check if kmer_index is valid and look up position
	if (kmer_index < this->current_index) {
		bool allele_id = this->path_to_allele.at(path_index);
		return (this->alleles.at(allele_id).kmer_path.get_position(kmer_index) > 0);
	} else {
		throw runtime_error("BiallelicUniqueKmers::kmer_on_path: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}


bool BiallelicUniqueKmers::kmer_on_allele(size_t kmer_index, size_t allele_id) const {
	return this->alleles.at(allele_id).kmer_path.get_position(kmer_index);
}


unsigned short BiallelicUniqueKmers::get_readcount_of(size_t kmer_index) {
	if (kmer_index < this->current_index) {
		return this->kmer_to_count[kmer_index];
	} else {
		throw runtime_error("BiallelicUniqueKmers::get_readcount_of: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

void BiallelicUniqueKmers::update_readcount(size_t kmer_index, unsigned short new_count) {
	if (kmer_index <  this->current_index) {
		this->kmer_to_count[kmer_index] = new_count;
	} else {
		throw runtime_error("BiallelicUniqueKmers::update_readcount: requested kmer index: " + to_string(kmer_index) + " does not exist.");
	}
}

size_t BiallelicUniqueKmers::size() const {
	return this->current_index;
}

unsigned short BiallelicUniqueKmers::get_nr_paths() const {
	return this->path_to_allele.size();
}

void BiallelicUniqueKmers::get_path_ids(vector<unsigned short>& p, vector<unsigned short>& a, vector<unsigned short>* only_include) {
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

void BiallelicUniqueKmers::get_allele_ids(vector<unsigned short>& a) {
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		a.push_back(it->first);
	}
}


void BiallelicUniqueKmers::get_defined_allele_ids(std::vector<unsigned short>& a) {
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		if (!it->second.is_undefined) a.push_back(it->first);
	}
}


ostream& operator<< (ostream& stream, const BiallelicUniqueKmers& uk) {
	stream << "BiallelicUniqueKmers for variant: " << uk.variant_pos << endl;
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


map<unsigned short, int> BiallelicUniqueKmers::kmers_on_alleles () const {
	map<unsigned short, int> result;
	for (auto it = this->alleles.begin(); it != this->alleles.end(); ++it) {
		result[it->first] = alleles.at(it->first).kmer_path.nr_kmers();
	}
	return result;
}


unsigned short BiallelicUniqueKmers::kmers_on_allele(unsigned short allele_id) const {
	if ((allele_id != 1) && (allele_id != 0) ) {
		throw runtime_error("BiallelicUniqueKmers::kmers_on_allele: allele_id must be either 0 or 1.");
	}
	return alleles.at(allele_id).kmer_path.nr_kmers();
}


unsigned short BiallelicUniqueKmers::present_kmers_on_allele(unsigned short allele_id) const {
	if ((allele_id != 1) && (allele_id != 0) ) {
		throw runtime_error("BiallelicUniqueKmers::present_kmers_on_allele: allele_id must be either 0 or 1.");
	}
	unsigned short result = 0;
	for (size_t i = 0; i < this->kmer_to_count.size(); ++i) {
		if (kmer_to_count[i] < 3) continue;
		if (alleles.at(allele_id).kmer_path.get_position(i) > 0) result += 1;
	}
	return result;
}

float BiallelicUniqueKmers::fraction_present_kmers_on_allele(unsigned short allele_id) const {
	if ((allele_id != 1) && (allele_id != 0) ) {
		throw runtime_error("BiallelicUniqueKmers::fraction_present_kmers_on_allele: allele_id must be either 0 or 1.");
	}
	unsigned short total = this->kmers_on_allele(allele_id);
	if (total > 0) return this->present_kmers_on_allele(allele_id) / (float) total;
	return 1.0;
}

bool BiallelicUniqueKmers::is_undefined_allele (unsigned short allele_id) const {
	if ((allele_id != 1) && (allele_id != 0) ) {
		throw runtime_error("BiallelicUniqueKmers::is_undefined_allele: allele_id must be either 0 or 1.");
	}
	// check if allele id exists
	auto it = this->alleles.find(allele_id);
	if (it != this->alleles.end()) {
		return it->second.is_undefined;
	} else {
		return false;
	}
}

void BiallelicUniqueKmers::set_undefined_allele (unsigned short allele_id) {
	if ((allele_id != 1) && (allele_id != 0) ) {
		throw runtime_error("BiallelicUniqueKmers::set_undefined_allele: allele_id " + to_string(allele_id) + " does not exist.");
	}
	auto it = this->alleles.find(allele_id);
	if (it == this->alleles.end()) {
		throw runtime_error("BiallelicUniqueKmers::set_undefined_allele: allele_id " + to_string(allele_id) + " does not exist.");
	}
	this->alleles[allele_id].is_undefined = true;
}

unsigned short BiallelicUniqueKmers::get_allele(unsigned short path_id) const {
	if (path_id >= this->path_to_allele.size()) {
		throw runtime_error("BiallelicUniqueKmers:get_allele: index out of bounds.");
	}
	return this->path_to_allele[path_id];
}


void BiallelicUniqueKmers::update_paths(vector<unsigned short>& path_ids) {
	size_t nr_paths = path_ids.size();
	vector<bool> updated_path_to_allele(nr_paths);
	map<bool, AlleleInfo> updated_alleles;
	vector<bool> undefined_alleles;

	for (size_t i = 0; i < path_ids.size(); ++i) {
		bool allele = this->path_to_allele[path_ids[i]];
		updated_path_to_allele[i] = allele;
		updated_alleles[allele] = this->alleles[allele];
	}

	// update the KmerPath objects and the kmer counts
	map<size_t, vector<unsigned short>> kmer_to_alleles;
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

void BiallelicUniqueKmers::print_kmer_matrix(string chromosome) const {
	for (auto a : this->alleles) {
		cout << chromosome << "\t" << this->variant_pos << "\t" << a.second.kmer_path << endl;
	}
}