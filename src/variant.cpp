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

Variant::Variant(string left_flank, string right_flank, string chromosome, size_t start_position, size_t end_position, vector<string> alleles, vector<unsigned char> paths)
	:left_flank(left_flank),
	 right_flank(right_flank),
	 chromosome(chromosome),
	 start_positions({start_position}),
	 end_positions({end_position}),
	 paths(paths),
	 flanks_added(false)
{
	for (auto a : alleles) {
		this->alleles.push_back({DnaSequence(a)});
	}
	this->set_values();
}

Variant::Variant(DnaSequence& left_flank, DnaSequence& right_flank, string chromosome, size_t start_position, size_t end_position, vector<DnaSequence>& alleles, vector<unsigned char>& paths)
	:left_flank(left_flank),
	 right_flank(right_flank),
	 chromosome(chromosome),
	 start_positions({start_position}),
	 end_positions({end_position}),
	 paths(paths),
	 flanks_added(false)
{
	for (auto a : alleles) {
		this->alleles.push_back({a});
	}
	this->set_values();
}

void Variant::set_values() {
	// find out which alleles are not covered by any paths
	vector<DnaSequence> uncovered;
	assert(alleles.size() < 256);
	for (unsigned char i = 0; i < this->alleles.size(); ++i) {
		if (find(this->paths.begin(), this->paths.end(), i) == this->paths.end()) {
			// allele not covered
			uncovered.push_back(this->alleles[i][0]);
		}
	}
	this->uncovered_alleles.push_back(uncovered);
	size_t start_position = this->start_positions[0];
	size_t end_position = this->end_positions[0];

	// check if flanks have same length
	if (this->left_flank.size() != this->right_flank.size()){
		throw runtime_error("Variant::Variant: left and right flanks have different sizes.");
	}

	// check if start and end positions are valid
	if (end_position <= start_position) {
		throw runtime_error("Variant::Variant: end position is smaller or equal to start position.");
	}

	// check if length of ref allele matches end position
	size_t ref_len = this->alleles[0][0].size();
	if (ref_len != (end_position-start_position)) {
		throw runtime_error("Variant::Variant: end position does not match length of reference allele.");
	}

	// check if paths are valid
	size_t nr_alleles = this->alleles.size();
	for (auto p : this->paths) {
		if (p >= nr_alleles) {
			throw runtime_error("Variant::Variant: allele ids given in paths are invalid.");
		}
	}
}

void Variant::add_flanking_sequence(){
	if (this->flanks_added) return;
	for (size_t i = 0; i < this->alleles.size(); ++i) {
		this->alleles[i].insert(this->alleles[i].begin(), this->left_flank);
		this->alleles[i].push_back(this->right_flank);
	}
	this->flanks_added = true;
}

void Variant::remove_flanking_sequence() {
	if (!this->flanks_added) return;
	for (size_t i = 0; i < this->alleles.size(); ++i) {
		vector<DnaSequence> allele = this->alleles[i];
		allele.erase(allele.begin());
		allele.erase(allele.end()-1);
		this->alleles[i] = allele;
	}
	this->flanks_added = false;
}

size_t Variant::nr_of_alleles() const {
	return this->alleles.size();
}

size_t Variant::nr_of_paths() const {
	return this->paths.size();
}

string Variant::get_allele_string(size_t index) const {
	if (index < this->alleles.size()) {
		DnaSequence result;
		for(auto s : this->alleles.at(index)) {
			result.append(s);
		}
		return result.to_string();
	} else {
		throw runtime_error("Variant::get_allele_string: Index out of bounds.");
	}
}

DnaSequence Variant::get_allele_sequence(size_t index) const {
	if (index < this->alleles.size()) {
		DnaSequence result;
		for(auto s : this->alleles.at(index)) {
			result.append(s);
		}
		return result;
	} else {
		throw runtime_error("Variant::get_allele_sequence: Index out of bounds.");
	}
}

size_t Variant::get_start_position() const {
	return this->start_positions[0];
}

size_t Variant::get_end_position() const {
	return this->end_positions[ this->end_positions.size() - 1];
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
	if (v2.get_start_position() < this->get_end_position()){
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
	unsigned int dist = v2.get_start_position() - this->get_end_position();
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
	vector<vector<DnaSequence>> new_alleles;
	unsigned char allele_index = 0;
	for (auto it = path_to_index.begin(); it != path_to_index.end(); ++it) {
		assert(allele_index < 256);
		for (auto e : it->second) {
			new_paths[e] = allele_index;
		}
		vector<DnaSequence> left_allele = this->alleles.at(it->first.first);
		vector<DnaSequence> right_allele = v2.alleles.at(it->first.second);
		DnaSequence flank;
		this->right_flank.substr(0, v2.get_start_position() - this->get_end_position(), flank);
		left_allele.push_back(flank);
		left_allele.insert(left_allele.end(), right_allele.begin(), right_allele.end());
		new_alleles.push_back(left_allele);
		allele_index += 1;
	}

	// update variant
	this->start_positions.insert(this->start_positions.end(), v2.start_positions.begin(), v2.start_positions.end());
	this->right_flank = v2.right_flank;
	this->end_positions.insert(this->end_positions.end(), v2.end_positions.begin(), v2.end_positions.end());
	this->alleles = new_alleles;
	this->uncovered_alleles.insert(this->uncovered_alleles.end(), v2.uncovered_alleles.begin(), v2.uncovered_alleles.end());
	this->paths = new_paths;
}

void Variant::separate_variants (vector<Variant>* resulting_variants, const GenotypingResult* input_genotyping, vector<GenotypingResult>* resulting_genotyping) const {
	size_t nr_variants = this->start_positions.size();
	assert (this->uncovered_alleles.size() == nr_variants);
	// collect the allele sequences for all variants 
	vector<vector<DnaSequence>> alleles_per_variant (nr_variants);
	for (size_t i = 0; i < this->alleles.size(); ++i) {
		vector<DnaSequence> allele = this->alleles.at(i);
		size_t start_index = 0;
		if (this->flanks_added) start_index = 1;
		for (size_t a = 0; a < nr_variants; a++) {
			size_t index = a*2 + start_index;
			DnaSequence sequence = allele.at(index);
			auto it = find(alleles_per_variant.at(a).begin(), alleles_per_variant.at(a).end(), sequence);
			if (it == alleles_per_variant.at(a).end()) {
				alleles_per_variant.at(a).push_back(sequence);
			}
		}
	}

	// construct paths
	vector<vector<unsigned char>> paths_per_variant (nr_variants);
	for (size_t i = 0; i < this->paths.size(); ++i) {
		unsigned char a = this->get_allele_on_path(i);
		vector<DnaSequence> allele = this->alleles.at(a);
		size_t start_index = 0;
		if (this->flanks_added) start_index = 1;
		for (size_t v = 0; v < nr_variants; v++) {
			// get allele at this position
			size_t index = v*2 + start_index;
			DnaSequence sequence = allele.at(index);
			// get index of this allele
			auto it = find(alleles_per_variant.at(v).begin(), alleles_per_variant.at(v).end(), sequence);
			unsigned char allele_id = distance(alleles_per_variant.at(v).begin(), it);
			paths_per_variant.at(v).push_back(allele_id);
		}
	}

	// use reference allele to construct flanking sequences for each variant
	vector<DnaSequence> reference_allele = this->alleles.at(0);
	if (!this->flanks_added) {
		reference_allele.insert(reference_allele.begin(), this->left_flank);
		reference_allele.push_back(this->right_flank);
	}

	for (size_t i = 0; i < nr_variants; ++i) {
		size_t index = i*2;
		if (this->flanks_added) index += 1;
		DnaSequence left = construct_left_flank(reference_allele, i*2 + 1, this->left_flank.size());
		DnaSequence right = construct_right_flank(reference_allele, i*2 + 1, this->right_flank.size());

		// construct variant
		vector<DnaSequence> uncovered = this->uncovered_alleles[i];
		// add uncovered alleles to list of alleles
		vector<DnaSequence> new_alleles = alleles_per_variant.at(i);
		for (auto allele : uncovered) {
			// only add allele if it is not yet in list of alleles
			if (find(new_alleles.begin(), new_alleles.end(), allele) == new_alleles.end()) {
				new_alleles.push_back(allele);
			}
		}
		// construct new variant object
		Variant v(left, right, this->chromosome, this->start_positions.at(i), this->end_positions.at(i), new_alleles, paths_per_variant.at(i));
		resulting_variants->push_back(v);

		if (input_genotyping != nullptr) {
			// construct GenotypingResult
			GenotypingResult g;
			// iterate through all genotypes and determine the genotype likelihoods for single variant
			for (size_t a0 = 0; a0 < this->nr_of_alleles(); ++a0) {
				for (size_t a1 = a0; a1 < this->nr_of_alleles(); ++a1) {
					// determine alleles genotype corresponds to
					vector<DnaSequence> allele0 = this->alleles.at(a0);
					vector<DnaSequence> allele1 = this->alleles.at(a1);
					// determine alleles for current variant
					auto it0 = find(alleles_per_variant.at(i).begin(), alleles_per_variant.at(i).end(), allele0.at(index));
					unsigned char single_allele0 = distance(alleles_per_variant.at(i).begin(), it0);
					auto it1 = find(alleles_per_variant.at(i).begin(), alleles_per_variant.at(i).end(), allele1.at(index));
					unsigned char single_allele1 = distance(alleles_per_variant.at(i).begin(), it1);
					// update genotype likelihood
					long double combined_likelihood = input_genotyping->get_genotype_likelihood(a0, a1);
					g.add_to_likelihood(single_allele0, single_allele1, combined_likelihood);
				}
			}
			// get the haplotype alleles of the combined variant
			pair<unsigned char,unsigned char> haplotype = input_genotyping->get_haplotype();
			// get corresponding alleles for current variant
			vector<DnaSequence> allele0 = this->alleles.at(haplotype.first);
			vector<DnaSequence> allele1 = this->alleles.at(haplotype.second);
			auto it0 = find(alleles_per_variant.at(i).begin(), alleles_per_variant.at(i).end(), allele0.at(index));
			unsigned char single_haplotype0 = distance(alleles_per_variant.at(i).begin(), it0);
			auto it1 = find(alleles_per_variant.at(i).begin(), alleles_per_variant.at(i).end(), allele1.at(index));
			unsigned char single_haplotype1 = distance(alleles_per_variant.at(i).begin(), it1);
			// update result
			g.add_first_haplotype_allele(single_haplotype0);
			g.add_second_haplotype_allele(single_haplotype1);
			g.set_nr_unique_kmers(input_genotyping->get_nr_unique_kmers());
			resulting_genotyping->push_back(g);
		}
	}
}

bool Variant::is_combined() const {
	return (this->start_positions.size() > 1);
}

ostream& operator<<(ostream& os, const Variant& var) {
	os << "left flank:\t" << var.left_flank.to_string() << endl;
	os << "right flank:\t" << var.right_flank.to_string() << endl;
	os << "position:\t" << var.chromosome << ":" << var.start_positions[0] << "-" << var.end_positions[0] << endl;
	os << "alleles:" << endl;
	for (size_t i = 0; i < var.alleles.size(); ++i) {
		os << i << ":\t";
		for (auto s : var.alleles[i]) {
			os << s.to_string();
		}
		os << endl;
	}
	os << "alleles not covered by any path:" << endl;
	for (size_t i = 0; i < var.uncovered_alleles.size(); ++i) {
		os << "{";
		for (size_t j = 0; j < var.uncovered_alleles[i].size(); ++j) {
			if (j > 0) os << ",";
			os << var.uncovered_alleles[i][j].to_string();
		}
		os << "}" << endl;
	}
	os << "paths:" << endl;
	for (size_t i = 0; i < var.paths.size(); ++i) {
		os << (size_t) var.paths[i] << "\t";
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
	if (v1.start_positions != v2.start_positions) return false;
	if (v1.end_positions != v2.end_positions) return false;

	// check alleles
	if (v1.alleles != v2.alleles) return false;

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
