#include "uniquekmercomputer.hpp"
#include <iostream>
#include <map>

using namespace std;

void unique_kmers(string allele, size_t index, size_t kmer_size, map<string, vector<size_t>>& occurences) {
	// enumerate kmers
	map <string, size_t> counts;
	for (size_t i = 0; i <= (allele.length() - kmer_size); ++i) {
		string kmer = allele.substr(i, kmer_size);
		counts[kmer] += 1;
	}
	// determine kmers unique to allele
	for (auto const& entry : counts) {
		if (entry.second == 1) occurences[entry.first].push_back(index);
	}
}

UniqueKmerComputer::UniqueKmerComputer (KmerCounter* genomic_kmers, KmerCounter* read_kmers, VariantReader* variants, string chromosome, size_t kmer_coverage)
	:genomic_kmers(genomic_kmers),
	 read_kmers(read_kmers),
	 variants(variants),
	 chromosome(chromosome),
	 kmer_coverage(kmer_coverage)
{
	double cn0 = 0.1;
	double cn1 = kmer_coverage / 2.0;
	double cn2 = kmer_coverage;
	cerr << "Using Poisson distributions with parameters: " << cn0 << " " << cn1 << " " << cn2 << endl;
	this->probability.set_parameters(cn0, cn1 , cn2);
}


// need to know which kmer is on which allele to be able to identify which kmer is on which path
vector<UniqueKmers>* UniqueKmerComputer::compute_unique_kmers() const {
	vector<UniqueKmers>* result = new vector<UniqueKmers>;
	size_t nr_variants = this->variants->size_of(this->chromosome);
	for (size_t v = 0; v < nr_variants; ++v) {
		map <string, vector<size_t>> occurences;
		const Variant& variant = this->variants->get_variant(this->chromosome, v);
		UniqueKmers u (v, variant.get_start_position());
		size_t nr_alleles = variant.nr_of_alleles();

		// insert empty paths (to also capture paths for which no unique kmers exist)
		for (size_t p = 0; p < variant.nr_of_paths(); ++p) {
			u.insert_empty_path(p);
		}
		for (size_t a = 0; a < nr_alleles; ++a) {
			// enumerate kmers and identify those with copynumber 1
			string allele = variant.get_allele_sequence(a);
			size_t kmer_size = this->variants->get_kmer_size();
			unique_kmers(allele, a, kmer_size, occurences);
		}

		// check if kmers occur elsewhere in the genome
		for (auto& kmer : occurences) {
			size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
			size_t local_count = kmer.second.size();

			if ( (genomic_count - local_count) == 0 ) {
				// kmer unique to this region
				// determine read kmercount for this kmer
				size_t read_kmercount = this->read_kmers->getKmerAbundance(kmer.first);

				// determine on which paths kmer occurs
				vector<size_t> paths;
				for (auto& allele : kmer.second) {
					variant.get_paths_of_allele(allele, paths);
				}

				// skip kmers with "too extreme" counts
				// TODO: value ok?
				if (read_kmercount > (this->kmer_coverage)) {
					continue;
				}

				// determine probabilities
				long double p_cn0 = this->probability.get_probability(0, read_kmercount);
				long double p_cn1 = this->probability.get_probability(1, read_kmercount);
				long double p_cn2 = this->probability.get_probability(2, read_kmercount);

				// skip kmers with only 0 probabilities
				if ( (p_cn0 > 0) || (p_cn1 > 0) || (p_cn2 > 0) ) {
					CopyNumber cn(p_cn0, p_cn1, p_cn2);
					// construct UniqueKmers object
					u.insert_kmer(cn, paths);
				}
			}
		}
		result->push_back(u);
	}
	return result;
}
