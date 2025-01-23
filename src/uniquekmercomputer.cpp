#include "uniquekmercomputer.hpp"
#include <iostream>
#include <cassert>
#include <map>
#include <queue>

using namespace std;

void unique_kmers(DnaSequence& allele, unsigned short index, size_t kmer_size, map<jellyfish::mer_dna, vector<unsigned short>>& occurences) {
	//enumerate kmers
	map<jellyfish::mer_dna, size_t> counts;
	size_t extra_shifts = kmer_size;
	jellyfish::mer_dna::k(kmer_size);
	jellyfish::mer_dna current_kmer("");
	for (size_t i = 0; i < allele.size(); ++i) {
		char current_base = allele[i];
		if (extra_shifts == 0) {
			counts[current_kmer] += 1;
		}
		if (  ( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
			extra_shifts = kmer_size + 1;
		}
		current_kmer.shift_left(current_base);
		if (extra_shifts > 0) extra_shifts -= 1;
	}
	counts[current_kmer] += 1;

	// determine kmers unique to allele
	for (auto const& entry : counts) {
		if (entry.second == 1) occurences[entry.first].push_back(index);
	}
}

UniqueKmerComputer::UniqueKmerComputer (KmerCounter* genomic_kmers, shared_ptr<KmerCounter> read_kmers, shared_ptr<Graph> variants, size_t kmer_coverage)
	:genomic_kmers(genomic_kmers),
	 read_kmers(read_kmers),
	 variants(variants),
	 chromosome(variants->get_chromosome()),
	 kmer_coverage(kmer_coverage)
{
	jellyfish::mer_dna::k(this->variants->get_kmer_size());
}


map<unsigned short, vector<jellyfish::mer_dna>> UniqueKmerComputer::select_kmers(const Variant* variant, std::map <jellyfish::mer_dna, vector<unsigned short>>& occurences) {

	size_t nr_selected = 0;
	map<unsigned short, vector<jellyfish::mer_dna>> result;
	map<unsigned short, queue<jellyfish::mer_dna>> allele_to_kmers;
	for (auto kmer : occurences){

		size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
		size_t local_count = kmer.second.size();

		// if kmer is not unique to the region, skip it
		if ( (genomic_count - local_count) != 0 ) continue;

		// if kmer occurs on more than one allele, skip it
		if (local_count > 1) continue;


		// skip alleles not covered by any paths
		assert(local_count == 1);
		vector<size_t> paths;
		variant->get_paths_of_allele(kmer.second[0], paths);
		if (paths.size() == 0) continue;

		allele_to_kmers[kmer.second[0]].push(kmer.first);
	}

	bool keep_adding = true;
	while ( (nr_selected < 301) && (keep_adding) ) {
		bool kmer_added = false;
		for (auto a : allele_to_kmers) {
			// pick at most 32 kmers per allele
			if ( (a.second.size() > 0) && (result[a.first].size() < 32)) {
				result[a.first].push_back(a.second.front());
				a.second.pop();
				kmer_added = true;
				nr_selected += 1;
			}
		}
		keep_adding = kmer_added;
	}
	return result;
}


void UniqueKmerComputer::compute_unique_kmers(vector<shared_ptr<UniqueKmers>>* result, ProbabilityTable* probabilities, bool delete_processed_variants) {
	size_t nr_variants = this->variants->size();
	for (size_t v = 0; v < nr_variants; ++v) {

		// set parameters of distributions
		size_t kmer_size = this->variants->get_kmer_size();
		double kmer_coverage = compute_local_coverage(v, 2*kmer_size);
		
		map <jellyfish::mer_dna, vector<unsigned short>> occurences;
		const Variant& variant = this->variants->get_variant(v);
	
		vector<unsigned short> path_to_alleles;
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned short a = variant.get_allele_on_path(p);
			path_to_alleles.push_back(a);
		}

		shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers>(new UniqueKmers(variant.get_start_position(), path_to_alleles));
		u->set_coverage(kmer_coverage);
		size_t nr_alleles = variant.nr_of_alleles();

		for (unsigned short a = 0; a < nr_alleles; ++a) {
			// consider all alleles not undefined
			if (variant.is_undefined_allele(a)) {
				// skip kmers of alleles that are undefined
				u->set_undefined_allele(a);
				continue;
			}
			DnaSequence allele = variant.get_allele_sequence(a);
			unique_kmers(allele, a, kmer_size, occurences);
		}

		// select unique kmers to be used
		map<unsigned short, vector<jellyfish::mer_dna>> allele_to_kmers = select_kmers(&variant, occurences);

		// construct UniqueKmers object
		for (auto& a : allele_to_kmers) {
			for (auto& kmer : a.second) {
				size_t read_kmercount = this->read_kmers->getKmerAbundance(kmer);
				CopyNumber cn = probabilities->get_probability(kmer_coverage, read_kmercount);
				long double p_cn0 = cn.get_probability_of(0);
				long double p_cn1 = cn.get_probability_of(1);
				long double p_cn2 = cn.get_probability_of(2);

				// skip kmers with only 0 probabilities
				if ( (p_cn0 > 0) || (p_cn1 > 0) || (p_cn2 > 0) ) {
					// insert the kmer
					vector<unsigned short> alleles = {a.first};
					u->insert_kmer(read_kmercount, alleles);
				}
			}
		}

		result->push_back(u);

		if (delete_processed_variants) {
			if (v > 0) {
				// previous variant object no longer needed
				this->variants->delete_variant(v - 1);
			}
			if (v == (nr_variants - 1)) {
				// last variant object, can be deleted
				this->variants->delete_variant(v);
			}
		}
	}
}

void UniqueKmerComputer::compute_empty(vector<UniqueKmers*>* result) const {
	size_t nr_variants = this->variants->size();
	for (size_t v = 0; v < nr_variants; ++v) {
		const Variant& variant = this->variants->get_variant(v);
		vector<unsigned short> path_to_alleles;
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned short a = variant.get_allele_on_path(p);
			path_to_alleles.push_back(a);
		}
		UniqueKmers* u = new UniqueKmers(variant.get_start_position(), path_to_alleles);
		result->push_back(u);
	}
}

unsigned short UniqueKmerComputer::compute_local_coverage(size_t var_index, size_t length) {
	DnaSequence left_overhang;
	DnaSequence right_overhang;
	size_t total_coverage = 0;
	size_t total_kmers = 0;

	this->variants->get_left_overhang(var_index, length, left_overhang);
	this->variants->get_right_overhang(var_index, length, right_overhang);

	size_t kmer_size = this->variants->get_kmer_size();
	map <jellyfish::mer_dna, vector<unsigned short>> occurences;
	unique_kmers(left_overhang, 0, kmer_size, occurences);
	unique_kmers(right_overhang, 1, kmer_size, occurences);

	for (auto& kmer : occurences) {
		size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
		if (genomic_count == 1) {
			size_t read_count = this->read_kmers->getKmerAbundance(kmer.first);
			// ignore too extreme counts
			if ( (read_count < (this->kmer_coverage/4)) || (read_count > (this->kmer_coverage*4)) ) continue;
			total_coverage += read_count;
			total_kmers += 1;
		}
	}
	// in case no unique kmers were found, use constant kmer coverage
	if ((total_kmers > 0) && (total_coverage > 0)){
		return total_coverage / total_kmers;
	} else {
		return this->kmer_coverage;
	}
}
