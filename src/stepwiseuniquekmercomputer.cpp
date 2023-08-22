#include "stepwiseuniquekmercomputer.hpp"
#include <jellyfish/mer_dna.hpp>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>

using namespace std;

void stepwise_unique_kmers(DnaSequence& allele, unsigned char index, size_t kmer_size, map<jellyfish::mer_dna, vector<unsigned char>>& occurences) {
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

StepwiseUniqueKmerComputer::StepwiseUniqueKmerComputer (KmerCounter* genomic_kmers, shared_ptr<Graph> variants)
	:genomic_kmers(genomic_kmers),
	 variants(variants),
	 chromosome(variants->get_chromosome())
{
	jellyfish::mer_dna::k(this->variants->get_kmer_size());
}

void StepwiseUniqueKmerComputer::compute_unique_kmers(vector<shared_ptr<UniqueKmers>>* result, string filename , bool delete_processed_variants) {
	gzFile outfile = gzopen(filename.c_str(), "wb");
	if (!outfile) {
		stringstream ss;
		ss << "UniqueKmerComputer::compute_unique_kmers: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
		throw runtime_error(ss.str());
	}

	// write header of output file
	string header = "#chromosome\tstart\tend\tunique_kmers\tunique_kmers_overhang\n";
	gzwrite(outfile, header.c_str(), header.length());
	size_t kmer_size = this->variants->get_kmer_size();
	size_t overhang_size = 2*kmer_size;

	size_t nr_variants = this->variants->size();
	for (size_t v = 0; v < nr_variants; ++v) {
		// set parameters of distributions
		size_t kmer_size = this->variants->get_kmer_size();
		
		map <jellyfish::mer_dna, vector<unsigned char>> occurences;
		const Variant& variant = this->variants->get_variant(v);
		stringstream outline;
		outline << variant.get_chromosome() << "\t" << variant.get_start_position() << "\t" << variant.get_end_position() << "\t";
	
		vector<unsigned char> path_to_alleles;
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned char a = variant.get_allele_on_path(p);
			path_to_alleles.push_back(a);
		}

		shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers>(new UniqueKmers(variant.get_start_position(), path_to_alleles));
		// set for 0 for now, since we do not know the kmer coverage yet
		u->set_coverage(0);
		size_t nr_alleles = variant.nr_of_alleles();

		for (unsigned char a = 0; a < nr_alleles; ++a) {
			// consider all alleles not undefined
			if (variant.is_undefined_allele(a)) {
				// skip kmers of alleles that are undefined
				u->set_undefined_allele(a);
				continue;
			}
			DnaSequence allele = variant.get_allele_sequence(a);
			stepwise_unique_kmers(allele, a, kmer_size, occurences);
		}

		// check if kmers occur elsewhere in the genome
		size_t nr_kmers_used = 0;
		bool not_first = false;
		for (auto& kmer : occurences) {
			if (nr_kmers_used > 300) break;

			size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
			size_t local_count = kmer.second.size();
			if ( (genomic_count - local_count) == 0 ) {
				// kmer unique to this region
				// determine on which paths kmer occurs
				vector<size_t> paths;
				for (auto& allele : kmer.second) {
					variant.get_paths_of_allele(allele, paths);
				}

				// skip kmer that does not occur on any path (uncovered allele)
				if (paths.size() == 0) {
					continue;
				}

				// skip kmer that occurs on all paths (they do not give any information about a genotype)
				if (paths.size() == variant.nr_of_paths()) {
					continue;
				}

				// set read kmer count to 0 for now, since we don't know it yet
				u->insert_kmer(0, kmer.second);
				if (not_first) outline << ",";
				outline << kmer.first;
				not_first = true;
				nr_kmers_used += 1;

			}
		}

		// in case no kmers were written, print "nan"
		if (!not_first) outline << "nan";

		// write unique kmers of left and right overhang to file
		vector<string> flanking_kmers;
		determine_unique_flanking_kmers(v, overhang_size, flanking_kmers);
		not_first = false;
		outline << "\t";
		for (auto& kmer : flanking_kmers) {
			if (not_first) outline << ",";
			outline << kmer;
			not_first = true;
		}
		if (!not_first) outline << "nan";
		outline << endl;
		gzwrite(outfile, outline.str().c_str(), outline.str().size());

		result->push_back(u);

		// if requested, delete variant objects once they are no longer needed
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
	gzclose(outfile);
}


void StepwiseUniqueKmerComputer::compute_unique_kmers_fasta(vector<shared_ptr<UniqueKmers>>* result, string filename , bool delete_processed_variants) {
	gzFile outfile = gzopen(filename.c_str(), "wb");
	if (!outfile) {
		stringstream ss;
		ss << "UniqueKmerComputer::compute_unique_kmers: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
		throw runtime_error(ss.str());
	}

	size_t kmer_size = this->variants->get_kmer_size();
	size_t overhang_size = 2*kmer_size;

	size_t nr_variants = this->variants->size();
	for (size_t v = 0; v < nr_variants; ++v) {
		// set parameters of distributions
		size_t kmer_size = this->variants->get_kmer_size();
		
		map <jellyfish::mer_dna, vector<unsigned char>> occurences;
		const Variant& variant = this->variants->get_variant(v);
		stringstream outline;
	
		vector<unsigned char> path_to_alleles;
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned char a = variant.get_allele_on_path(p);
			path_to_alleles.push_back(a);
		}

		shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers>(new UniqueKmers(variant.get_start_position(), path_to_alleles));
		// set for 0 for now, since we do not know the kmer coverage yet
		u->set_coverage(0);
		size_t nr_alleles = variant.nr_of_alleles();

		for (unsigned char a = 0; a < nr_alleles; ++a) {
			// consider all alleles not undefined
			if (variant.is_undefined_allele(a)) {
				// skip kmers of alleles that are undefined
				u->set_undefined_allele(a);
				continue;
			}
			DnaSequence allele = variant.get_allele_sequence(a);
			stepwise_unique_kmers(allele, a, kmer_size, occurences);
		}

		// check if kmers occur elsewhere in the genome
		size_t nr_kmers_used = 0;
		for (auto& kmer : occurences) {
			if (nr_kmers_used > 300) break;

			size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
			size_t local_count = kmer.second.size();
			if ( (genomic_count - local_count) == 0 ) {
				// kmer unique to this region
				// determine on which paths kmer occurs
				vector<size_t> paths;
				for (auto& allele : kmer.second) {
					variant.get_paths_of_allele(allele, paths);
				}

				// skip kmer that does not occur on any path (uncovered allele)
				if (paths.size() == 0) {
					continue;
				}

				// skip kmer that occurs on all paths (they do not give any information about a genotype)
				if (paths.size() == variant.nr_of_paths()) {
					continue;
				}

				// set read kmer count to 0 for now, since we don't know it yet
				u->insert_kmer(0, kmer.second); // ><u/f>_<kmer-index>_<variant-index>_<variant-position>
				outline << ">u_" << nr_kmers_used << "_" << v << "_" << variant.get_start_position() << endl;
				outline << kmer.first << "\n";
				nr_kmers_used += 1;
			}
		}

		// write unique kmers of left and right overhang to file
		vector<string> flanking_kmers;
		determine_unique_flanking_kmers(v, overhang_size, flanking_kmers);
		size_t f_index = 0;
		for (auto& kmer : flanking_kmers) {
			outline << ">f_" << f_index << "_" << v << "_" << variant.get_start_position() << "\n";
			outline << kmer << "\n";
			f_index +=1;
		}
	
		gzwrite(outfile, outline.str().c_str(), outline.str().size());

		result->push_back(u);

		// if requested, delete variant objects once they are no longer needed
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
	gzclose(outfile);
}


void StepwiseUniqueKmerComputer::compute_empty(vector<shared_ptr<UniqueKmers>>* result) const {
	size_t nr_variants = this->variants->size();
	for (size_t v = 0; v < nr_variants; ++v) {
		const Variant& variant = this->variants->get_variant(v);
		vector<unsigned char> path_to_alleles;
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned char a = variant.get_allele_on_path(p);
			path_to_alleles.push_back(a);
		}
		shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers>(new UniqueKmers(variant.get_start_position(), path_to_alleles));
		result->push_back(u);
	}
}


void StepwiseUniqueKmerComputer::determine_unique_flanking_kmers(size_t var_index, size_t length, vector<string>& result) {
	DnaSequence left_overhang;
	DnaSequence right_overhang;

	this->variants->get_left_overhang(var_index, length, left_overhang);
	this->variants->get_right_overhang(var_index, length, right_overhang);

	size_t kmer_size = this->variants->get_kmer_size();
	map <jellyfish::mer_dna, vector<unsigned char>> occurences;
	stepwise_unique_kmers(left_overhang, 0, kmer_size, occurences);
	stepwise_unique_kmers(right_overhang, 1, kmer_size, occurences);

	for (auto& kmer : occurences) {
		size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer.first);
		if (genomic_count == 1) {
			result.push_back(kmer.first.to_str());
		}
	}
}