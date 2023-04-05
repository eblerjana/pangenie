#include <stdlib.h>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <map>
#include <cassert>
#include "variant.hpp"

using namespace std;

DnaSequence construct_left_flank(vector<DnaSequence>& alleles, size_t position, size_t length) {
	if ((alleles.size()-1) < position) {
		throw runtime_error("Variant::construct_left_flank: position too large.");
	}
	DnaSequence flank;
	for (int i = position-1; i >= 0; --i) {
		DnaSequence sequence = alleles.at(i);
		for (int j = sequence.size()-1; j >= 0; --j) {
			flank.append(sequence.base_at(j));
			if (flank.size() == length) break;
		}
		if (flank.size() == length) {
			break;
		}
	}
	flank.reverse();
	return flank;
}

DnaSequence construct_right_flank(vector<DnaSequence>& alleles, size_t position, size_t length) {
	if ((alleles.size()-1) < position) {
		throw runtime_error("Variant::construct_right_flank: position too large.");
	}
	DnaSequence flank;
	for (size_t i = position+1; i < alleles.size(); ++i) {
		DnaSequence sequence = alleles.at(i);
		for (size_t j = 0; j < sequence.size(); ++j) {
			flank.append(sequence.base_at(j));
			if (flank.size() == length) break;
		}
		if (flank.size() == length) {
			break;
		}
	}

	if (flank.size() != length) {
		throw runtime_error("Variant::construct_right_flank: not enough bases given at right side.");
	}
	return flank;
}

Variant::Variant(string left_flank, string right_flank, string chromosome, size_t start_position, size_t end_position, vector<string> alleles, vector<unsigned char> paths) //, string variant_id)
	:left_flank(left_flank),
	 right_flank(right_flank),
	 chromosome(chromosome),
	 start_position(start_position),
//	 variant_ids({variant_id}),
	 paths(paths),
	 flanks_added(false)

{
	if (alleles.size() > 255) {
		throw runtime_error("Variant::Variant: number of alleles per variant exceeds 256. Current implementation does not support higher numbers.");
	}

	if (paths.size() > 255) {
		throw runtime_error("Variant::Variant: number of paths exceeds 256. Current implementation does not support higher numbers.");
	}

	this->allele_sequences.push_back(vector<DnaSequence>());
	for (unsigned char i = 0; i < alleles.size(); ++i) {
		this->allele_sequences[0].push_back(DnaSequence(alleles[i]));
		this->allele_combinations.push_back({i});
	}

	this->set_values(end_position);
}

Variant::Variant(DnaSequence& left_flank, DnaSequence& right_flank, string chromosome, size_t start_position, size_t end_position, vector<DnaSequence>& alleles, vector<unsigned char>& paths) //, string variant_id)
	:left_flank(left_flank),
	 right_flank(right_flank),
	 chromosome(chromosome),
	 start_position(start_position),
//	 variant_ids({variant_id}),
	 paths(paths),
	 flanks_added(false)
{
	if (alleles.size() > 255) {
		throw runtime_error("Variant::Variant: number of alleles per variant exceeds 256. Current implementation does not support higher numbers.");
	}

	if (paths.size() > 255) {
		throw runtime_error("Variant::Variant: number of paths exceeds 256. Current implementation does not support higher numbers.");
	}

	this->allele_sequences.push_back(vector<DnaSequence>());
	for (unsigned char i = 0; i < alleles.size(); ++i) {
		this->allele_sequences[0].push_back(alleles[i]);
		this->allele_combinations.push_back({i});
	}

	this->set_values(end_position);
}

void Variant::set_values(size_t end_position) {
	// find out which alleles are not covered by any paths
	vector<unsigned char> uncovered;
	assert(allele_sequences.size() < 256);
	for (unsigned char i = 0; i < this->allele_sequences[0].size(); ++i) {
		if (find(this->paths.begin(), this->paths.end(), i) == this->paths.end()) {
			// allele not covered
			uncovered.push_back(i);
		}
	}

	this->uncovered_alleles.push_back(uncovered);

	// check if flanks have same length
	if (this->left_flank.size() != this->right_flank.size()){
		throw runtime_error("Variant::Variant: left and right flanks have different sizes.");
	}

	// check if start and end positions are valid
	if (end_position <= this->start_position) {
		throw runtime_error("Variant::Variant: end position is smaller or equal to start position.");
	}

	// check if length of ref allele matches end position
	size_t ref_len = this->allele_sequences[0][0].size();
	if (ref_len != (end_position - this->start_position)) {
		throw runtime_error("Variant::Variant: end position does not match length of reference allele.");
	}

	// check if paths are valid
	size_t nr_alleles = this->allele_sequences[0].size();
	for (auto p : this->paths) {
		if (p >= nr_alleles) {
			throw runtime_error("Variant::Variant: allele ids given in paths are invalid.");
		}
	}
}

void Variant::add_flanking_sequence(){
	this->flanks_added = true;
}

void Variant::remove_flanking_sequence() {
	this->flanks_added = false;
}

size_t Variant::nr_of_alleles() const {
	return this->allele_combinations.size();
}

size_t Variant::nr_of_paths() const {
	return this->paths.size();
}

string Variant::get_allele_string(size_t index) const {
	if (index < this->allele_combinations.size()) {
		DnaSequence result;
		if (this->flanks_added) {
			result = this->left_flank;
		}
		size_t nr_alleles = this->allele_combinations.at(index).size();
		for (size_t i = 0; i < nr_alleles; ++i) {
			result.append(this->allele_sequences[i][this->allele_combinations[index][i]]);
			if (i < (nr_alleles - 1)) {
				result.append(this->inner_flanks[i]);
			}
		}
		if (this->flanks_added) {
			result.append(this->right_flank);
		}
		return result.to_string();
	} else {
		throw runtime_error("Variant::get_allele_string: Index out of bounds.");
	}
}

DnaSequence Variant::get_allele_sequence(size_t index) const {
	if (index < this->allele_combinations.size()) {
		DnaSequence result;
		if (this->flanks_added) {
			result = this->left_flank;
		}
		size_t nr_alleles = this->allele_combinations.at(index).size();
		for (size_t i = 0; i < nr_alleles; ++i) {
			result.append(this->allele_sequences[i][this->allele_combinations[index][i]]);
			if (i < (nr_alleles - 1)) {
				result.append(this->inner_flanks[i]);
			}
		}
		if (this->flanks_added) {
			result.append(this->right_flank);
		}
		return result;
	} else {
		throw runtime_error("Variant::get_allele_string: Index out of bounds.");
	}
}

size_t Variant::get_start_position() const {
	return this->start_position;
}

size_t Variant::get_end_position() const {
	size_t end_position = this->start_position;
	for (size_t i = 0; i < this->allele_sequences.size(); ++i) {
		end_position += this->allele_sequences.at(i).at(0).size();
		if (i < (this->allele_sequences.size()-1)) {
			end_position += this->inner_flanks.at(i).size();
		}
	}
	return end_position;
}

string Variant::get_chromosome() const {
	return this->chromosome;
}

bool Variant::allele_on_path(unsigned char allele_index, size_t path_index) const {
	return (this->paths[path_index] == allele_index);
}

unsigned char Variant::get_allele_on_path(size_t path_index) const {
	return (this->paths[path_index]);
}

void Variant::get_paths_of_allele(unsigned char allele_index, std::vector<size_t>& result) const {
	for (size_t i = 0; i < this->paths.size(); ++i) {
		if (allele_on_path(allele_index, i)) {
			result.push_back(i);
		}
	}
}

void Variant::combine_variants (Variant const &v2){
	size_t end_position = this->get_end_position();
	if (v2.get_start_position() < end_position){
		throw runtime_error("Variant::combine_variants: Variants are overlapping.");
	}
	if (this->flanks_added || v2.flanks_added){
		throw runtime_error("Variant::combine_variants: Variant objects can only be combined if no flanks where added.");
	}
	size_t kmersize_v1 = this->left_flank.size();
	size_t kmersize_v2 = v2.left_flank.size();
	if (kmersize_v1 != kmersize_v2) {
		throw runtime_error("Variant::combine_variants: kmersizes are not the same.");
	}
	unsigned int dist = v2.get_start_position() - end_position;
	if ( (dist > kmersize_v1) || (this->chromosome != v2.chromosome) ){
		throw runtime_error("Variant::combine_variants: Variant objects are more that kmersize bases abart.");
	}
	if (this->paths.size() != v2.paths.size()){
		throw runtime_error("Variant::combine_variants: Variant objects not covered by the same paths.");
	}

	// consider all combinations of alleles defined by paths (=new alleles)
	map<size_t, pair<unsigned char,unsigned char>> index_to_path;
	map<pair<unsigned char,unsigned char>, vector<size_t>> path_to_index;

	for (size_t p = 0; p < this->paths.size(); ++p) {
		unsigned char left_allele = this->paths.at(p);
		unsigned char right_allele = v2.paths.at(p);
		index_to_path[p] = make_pair(left_allele, right_allele);
		path_to_index[make_pair(left_allele,right_allele)].push_back(p);
	}

	// add REF-REF allele
	pair<unsigned char,unsigned char> ref_path = make_pair(0,0);
	if (path_to_index.find(ref_path) == path_to_index.end()) {
		path_to_index[ref_path] = {};
	}
	vector<unsigned char> new_paths(this->paths.size());
	vector<vector<unsigned char>> new_alleles;
	unsigned char allele_index = 0;
	
	assert (path_to_index.size() < 256);
	
	// construct new allele sequences
	for (auto it = path_to_index.begin(); it != path_to_index.end(); ++it) {	
		for (auto e : it->second) {
			new_paths[e] = allele_index;
		}
		vector<unsigned char> left_allele = this->allele_combinations.at(it->first.first);
		vector<unsigned char> right_allele = v2.allele_combinations.at(it->first.second);
		left_allele.insert(left_allele.end(), right_allele.begin(), right_allele.end());
		new_alleles.push_back(left_allele);
		allele_index += 1;
	}

	// construct sequence between variants
	DnaSequence flank;
	this->right_flank.substr(0, v2.get_start_position() - end_position, flank);
	this->inner_flanks.push_back(flank);
	this->inner_flanks.insert(this->inner_flanks.end(), v2.inner_flanks.begin(), v2.inner_flanks.end());

	// update variant
	this->right_flank = v2.right_flank;
	this->allele_combinations = new_alleles;
	this->allele_sequences.insert(this->allele_sequences.end(), v2.allele_sequences.begin(), v2.allele_sequences.end());
	this->uncovered_alleles.insert(this->uncovered_alleles.end(), v2.uncovered_alleles.begin(), v2.uncovered_alleles.end());
	this->paths = new_paths;
//	this->variant_ids.insert(this->variant_ids.end(), v2.variant_ids.begin(), v2.variant_ids.end());
}

void Variant::separate_variants (vector<Variant>* resulting_variants, const GenotypingResult* input_genotyping, vector<GenotypingResult>* resulting_genotyping) const {
	size_t nr_variants = this->allele_sequences.size();
	assert (this->uncovered_alleles.size() == nr_variants);

	// construct paths
	vector<vector<unsigned char>> paths_per_variant (nr_variants);
	for (size_t i = 0; i < this->paths.size(); ++i) {
		unsigned char a = this->get_allele_on_path(i);
		vector<unsigned char> allele = this->allele_combinations.at(a);
		assert (allele.size() == nr_variants);
		for (size_t v = 0; v < nr_variants; v++) {
			// get allele at this position
			unsigned char allele_id = allele.at(v);
			paths_per_variant.at(v).push_back(allele_id);
		}
	}

	// use reference allele to construct flanking sequences for each variant
	vector<DnaSequence> reference_allele;
	for (size_t i = 0; i < nr_variants; ++i) {
		unsigned char allele_id = this->allele_combinations.at(0).at(i);
		reference_allele.push_back(this->allele_sequences.at(i).at(allele_id));
		if (i < (nr_variants - 1)) {
			reference_allele.push_back(this->inner_flanks[i]);
		}
	}

	reference_allele.insert(reference_allele.begin(), this->left_flank);
	reference_allele.push_back(this->right_flank);

	size_t current_start = this->start_position;

	for (size_t i = 0; i < nr_variants; ++i) {
		DnaSequence left = construct_left_flank(reference_allele, i*2 + 1, this->left_flank.size());
		DnaSequence right = construct_right_flank(reference_allele, i*2 + 1, this->right_flank.size());
		vector<DnaSequence> alleles = this->allele_sequences.at(i);
		size_t current_end = current_start + alleles[0].size();

		// construct new variant object
		Variant v(left, right, this->chromosome, current_start, current_end, alleles, paths_per_variant.at(i)); //, this->variant_ids.at(i));

		resulting_variants->push_back(v);
		if (input_genotyping != nullptr) {
			// construct GenotypingResult
			GenotypingResult g;
			// precompute alleles
			vector<unsigned char> precomputed_ids (this->nr_of_alleles());
			for (size_t a0 = 0; a0 < this->nr_of_alleles(); ++a0) {
				unsigned char single_allele0 = this->allele_combinations[a0][i];
				precomputed_ids[a0] = single_allele0;
			}
			// iterate through all genotypes and determine the genotype likelihoods for single variant
			for (size_t a0 = 0; a0 < this->nr_of_alleles(); ++a0) {
				// determine allele a0 genotype corresponds to
				unsigned char single_allele0 = precomputed_ids[a0];
				for (size_t a1 = a0; a1 < this->nr_of_alleles(); ++a1) {
					// determine allele a1 genotype corresponds to
					unsigned char single_allele1 = precomputed_ids[a1];
					// update genotype likelihood
					long double combined_likelihood = input_genotyping->get_genotype_likelihood(a0, a1);
					g.add_to_likelihood(single_allele0, single_allele1, combined_likelihood);
				}
			}
			// get the haplotype alleles of the combined variant
			pair<unsigned char,unsigned char> haplotype = input_genotyping->get_haplotype();
			// get corresponding alleles for current variant
			unsigned char single_haplotype0 = precomputed_ids[haplotype.first];
			unsigned char single_haplotype1 = precomputed_ids[haplotype.second];
			// update result
			g.add_first_haplotype_allele(single_haplotype0);
			g.add_second_haplotype_allele(single_haplotype1);
			resulting_genotyping->push_back(g);
		}
		// update start position
		current_start = current_end;
		if (i < (nr_variants-1)) {
			current_start += this->inner_flanks[i].size();
		}
	}
}


void Variant::variant_statistics (shared_ptr<UniqueKmers> unique_kmers, vector<VariantStats>& result) const {
	size_t nr_variants = this->allele_sequences.size();
	assert (this->uncovered_alleles.size() == nr_variants);

	for (size_t i = 0; i < nr_variants; ++i) {
		VariantStats v;
		// new allele -> unique kmer counts map
		map<unsigned char, int> new_kmer_counts;
		vector<unsigned char> precomputed_ids (this->nr_of_alleles());
		for (size_t a0 = 0; a0 < this->nr_of_alleles(); ++a0) {
			unsigned char single_allele0 = this->allele_combinations[a0][i];
			precomputed_ids[a0] = single_allele0;
		}
		// iterate through all alleles and determine number of unique kmers
		for (size_t a0 = 0; a0 < this->nr_of_alleles(); ++a0) {
			unsigned char single_allele0 = precomputed_ids[a0];
			// update unique kmer counts
			new_kmer_counts[single_allele0] += unique_kmers->kmers_on_alleles()[a0];
		}

		v.nr_unique_kmers = unique_kmers->size();
		v.coverage = unique_kmers->get_coverage();

		// determine ids of uncovered paths
		vector<unsigned char> uncovered_ids = this->uncovered_alleles[i];
		// set kmer counts of uncovered alleles to -1
		for (auto u: uncovered_ids) {
			new_kmer_counts[u] = -1;
		}
		
		v.kmer_counts = new_kmer_counts;
		result.push_back(v);
	}
}


bool Variant::is_combined() const {
	return (this->allele_sequences.size() > 1);
}

ostream& operator<<(ostream& os, const Variant& var) {
	os << "left flank:\t" << var.left_flank.to_string() << endl;
	os << "right flank:\t" << var.right_flank.to_string() << endl;	
	os << "position:\t" << var.chromosome << ":" << var.start_position << "-" << var.get_end_position() << endl;
	os << "alleles:" << endl;
	for (size_t i = 0; i < var.allele_combinations.size(); ++i) {
		os << i << ":\t";
		os << var.get_allele_string(i);
		os << endl;
	}
	os << "alleles not covered by any path:" << endl;
	for (size_t i = 0; i < var.uncovered_alleles.size(); ++i) {
		os << "{";
		for (size_t j = 0; j < var.uncovered_alleles[i].size(); ++j) {
			if (j > 0) os << ",";
			unsigned char id = var.uncovered_alleles[i][j];
			os << var.allele_sequences[i][id].to_string();
		}
		os << "}" << endl;
	}
	os << "paths:" << endl;
	for (size_t i = 0; i < var.paths.size(); ++i) {
		os << (size_t) var.paths[i] << "\t";
	}
	
	os << "inner flanks:" << endl;
	for (auto s : var.inner_flanks) {
		os << s.to_string() << endl;
	}
	return os;
}

bool operator==(const Variant& v1, const Variant& v2) {
	// check flanks
	if (v1.left_flank != v2.left_flank) return false;
	if (v1.right_flank != v2.right_flank) return false;

	// check chromosome
	if (v1.chromosome != v2.chromosome) return false;

	// check positions
	if (v1.start_position != v2.start_position) return false;
	if (v1.get_end_position() != v2.get_end_position()) return false;

	// check alleles
	if (v1.allele_sequences != v2.allele_sequences) return false;
	if (v1.allele_combinations != v2.allele_combinations) return false;
	if (v1.inner_flanks != v2.inner_flanks) return false;

	// check uncovered alleles
	if (v1.uncovered_alleles != v2.uncovered_alleles) return false;

	// check paths
	if (v1.paths != v2.paths) return false;

	// check flanks_added
	if (v1.flanks_added != v2.flanks_added) return false;

	return true;
}

bool operator!=(const Variant& v1, const Variant& v2) {
	return !(v1 == v2);
}

float Variant::allele_frequency(unsigned char allele_index, bool ignore_ref_path) const {
	if (this->paths.size() == 0) {
		return 0.0;
	}
	float freq = 0.0;
	for (auto a : this->paths) {
		if (a == allele_index) freq += 1;
	}
	unsigned int size = paths.size();
	if (ignore_ref_path) size -= 1.0;
	if (ignore_ref_path && (allele_index == 0)) {
		assert(freq >= 1.0);
		freq -= 1.0;
	}
	return freq / size;
}

string Variant::get_id() const {
//	string result = "";
//	for (size_t i = 0; i < this->variant_ids.size(); ++i) {
//		if (i > 0) result += ";";
//		result += this->variant_ids.at(i);
//	}
//	return result;
	return ".";
}

bool Variant::is_undefined_allele(size_t allele_id) const {
	DnaSequence allele_without_flanks;
	for (size_t i = 0; i < this->allele_combinations.at(allele_id).size(); ++i) {
		unsigned char allele = this->allele_combinations.at(allele_id).at(i);
		allele_without_flanks.append(this->allele_sequences.at(i).at(allele));
	}
	return allele_without_flanks.contains_undefined();
}

size_t Variant::nr_missing_alleles() const {
	size_t missing = 0;
	for (auto path : this->paths) {
		DnaSequence allele = this->get_allele_sequence(path);
		if (allele.contains_undefined()) missing += 1;
	}
	return missing;
}
