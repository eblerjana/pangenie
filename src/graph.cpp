#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <regex>
#include "graphbuilder.hpp"
#include "graph.hpp"

using namespace std;


string graph_get_date() {
	time_t t = time(0);
	tm* now = localtime(&t);
	ostringstream oss;
	oss << (now->tm_year+1900) << setfill('0') << setw(2) << (now->tm_mon+1) << setw(2) << now->tm_mday;
	return string(oss.str());
}

void Graph::insert_ids(vector<DnaSequence>& alleles, vector<string>& variant_ids, bool reference_added) {
	vector<unsigned char> index = graph_construct_index(alleles, reference_added);
	assert(index.size() < 256);
	// insert IDs in the lex. order of their corresponding alleles
	vector<string> sorted_ids;
	for (auto id : index) {
		sorted_ids.push_back(variant_ids[id]);
	}
	this->variant_ids.push_back(sorted_ids);
}

string Graph::get_ids(vector<string>& alleles, size_t variant_index, bool reference_added) {
	vector<unsigned char> index = graph_construct_index(alleles, reference_added);
	assert(index.size() < 256);
	vector<string> sorted_ids(index.size());
	for (unsigned char i = 0; i < index.size(); ++i) {
		sorted_ids[index[i]] = this->variant_ids.at(variant_index)[i];
	}
	string result = "";
	for (unsigned char i = 0; i < sorted_ids.size(); ++i) {
		if (i > 0) result += ',';
		result += sorted_ids[i];
	}
	return result;
}

Graph::Graph(FastaReader fasta_reader, string chromosome, size_t kmer_size, bool add_reference)
	: fasta_reader(fasta_reader),
	  chromosome(chromosome),
	  kmer_size(kmer_size),
	  add_reference(add_reference),
	  variants_deleted(false)
{}

size_t Graph::get_kmer_size() const {
	return this->kmer_size;
}

string Graph::get_chromosome() const {
	return this->chromosome;
}

size_t Graph::size() const {
	return this->variants.size();
}

void Graph::add_variant_cluster(vector<shared_ptr<Variant>>* cluster, vector<vector<string>>& variant_ids, bool only_defined_ids) {
	if (!cluster->empty()) {
		assert(cluster->size() == variant_ids.size());

		// store variant ids
		for (size_t i = 0; i < variant_ids.size(); ++i) {
			if (!variant_ids[i].empty()) {
				assert (cluster->at(i)->allele_sequences.size() == 1);
				if (only_defined_ids) {
					vector<DnaSequence> defined_alleles;
					for (auto a : cluster->at(i)->allele_sequences.at(0)) {
						if (a.contains_undefined()) continue;
						defined_alleles.push_back(a);
					}
					assert (defined_alleles.size() == (variant_ids[i].size()+1));
					insert_ids(defined_alleles, variant_ids[i], true);
				} else {
					insert_ids(cluster->at(i)->allele_sequences.at(0), variant_ids[i], true);
				}
			} else {
				this->variant_ids.push_back(vector<string>());
			}
		}

		// merge all variants in cluster
		shared_ptr<Variant> combined = cluster->at(0);
		for (size_t v = 1; v < cluster->size(); ++v) {
			combined->combine_variants(*cluster->at(v));
			cluster->operator[](v) = nullptr;
		}
		combined->add_flanking_sequence();
		this->variants.push_back(combined);

	}
}

const Variant& Graph::get_variant(size_t index) const {
	if (index < this->size()) {
		if (this->variants.at(index) != nullptr) {
			return *this->variants.at(index);
		} else {
			throw runtime_error("Graph::get_variant: variant was previously destroyed by delete_variant function.");
		}
	} else {
		throw runtime_error("Graph::get_variant: index out of bounds.");
	}
}

const FastaReader& Graph::get_fasta_reader() const {
	return this->fasta_reader;
}

void Graph::write_genotypes(string filename, const vector<GenotypingResult>& genotyping_result, bool write_header, string sample, bool ignore_imputed) {
	if (this->variants_deleted) {
		throw runtime_error("Graph::write_genotypes_of: variants have been deleted by delete_variant funtion. Re-build object.");
	}

	if (genotyping_result.size() != this->size()) {
		throw runtime_error("Graph::write_genotypes_of: number of variants and number of computed genotypes differ.");
	}

	ofstream genotyping_outfile;
	if (write_header) {
		genotyping_outfile.open(filename);
		if (! genotyping_outfile.is_open()) {
			throw runtime_error("Graph::write_genotypes_of: genotyping output file cannot be opened. Note that the filename must not contain non-existing directories.");
		}
		
		// write VCF header lines
		genotyping_outfile << "##fileformat=VCFv4.2" << endl;
		genotyping_outfile << "##fileDate=" << graph_get_date() << endl;
		// TODO output command line
		genotyping_outfile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
		genotyping_outfile << "##INFO=<ID=UK,Number=1,Type=Integer,Description=\"Total number of unique kmers.\">" << endl;
		genotyping_outfile << "##INFO=<ID=AK,Number=R,Type=Integer,Description=\"Number of unique kmers per allele. Will be -1 for alleles not covered by any input haplotype path\">" << endl;
		genotyping_outfile << "##INFO=<ID=MA,Number=1,Type=Integer,Description=\"Number of alleles missing in panel haplotypes.\">" << endl;
		genotyping_outfile << "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">" << endl;
		genotyping_outfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
		genotyping_outfile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality: phred scaled probability that the genotype is wrong.\">" << endl;
		genotyping_outfile << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Comma-separated log10-scaled genotype likelihoods for absent, heterozygous, homozygous.\">" << endl;
		genotyping_outfile << "##FORMAT=<ID=KC,Number=1,Type=Float,Description=\"Local kmer coverage.\">" << endl;
		genotyping_outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample << endl; 
	} else {
		genotyping_outfile.open(filename, std::ios_base::app);
		if (! genotyping_outfile.is_open()) {
			throw runtime_error("Graph::write_genotypes_of: genotyping output file cannot be opened. Note that the filename must not contain non-existing directories.");
		}
	}

	size_t counter = 0;
	for (size_t i = 0; i < this->size(); ++i) {
		shared_ptr<Variant> variant = this->variants.at(i);
		unsigned short coverage = genotyping_result.at(i).coverage();
		unsigned short nr_unique_kmers = genotyping_result.at(i).nr_unique_kmers();

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
		variant->separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);

		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			genotyping_outfile << v.get_chromosome() << "\t"; // CHROM
			genotyping_outfile << (v.get_start_position() + 1) << "\t"; // POS
			genotyping_outfile << v.get_id() << "\t"; // ID
			genotyping_outfile << v.get_allele_string(0) << "\t"; // REF

			// get alternative allele
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "Graph::write_genotypes_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			vector<string> alt_alleles;
			vector<unsigned char> defined_alleles = {0};
			for (size_t i = 1; i < nr_alleles; ++i) {
				DnaSequence allele = v.get_allele_sequence(i);
				// skip alleles that are undefined
				if (!v.is_undefined_allele(i)) {
					alt_alleles.push_back(allele.to_string());
					defined_alleles.push_back(i);
				}
			}

			string alt_string = "";
			for (unsigned char a = 0; a < alt_alleles.size(); ++a) {
				if (a > 0) alt_string += ',';
				alt_string += alt_alleles[a];
			}

			genotyping_outfile << alt_string << "\t"; // ALT
			genotyping_outfile << ".\t"; // QUAL
			genotyping_outfile << "PASS" << "\t"; // FILTER
			// output allele frequencies of all alleles
			ostringstream info;
			info << "AF="; // AF
			for (unsigned int a = 1; a < defined_alleles.size(); ++a) {
				if (a > 1) info << ",";
				info << setprecision(6) << v.allele_frequency(defined_alleles[a], this->add_reference);				
			}

			// keep only likelihoods for genotypes with defined alleles
			size_t nr_missing = v.nr_missing_alleles();
			GenotypingResult genotype_likelihoods = singleton_likelihoods.at(j);
			if (nr_missing > 0) genotype_likelihoods = singleton_likelihoods.at(j).get_specific_likelihoods(defined_alleles);
			nr_alleles = defined_alleles.size();

			info << ";UK=" << nr_unique_kmers; // UK
			info << ";MA=" << nr_missing;
	
			// if IDs were given in input, write them to output as well
			if (!this->variant_ids[counter].empty()) info << ";ID=" << get_ids(alt_alleles, counter, false);
			genotyping_outfile << info.str() << "\t"; // INFO
			genotyping_outfile << "GT:GQ:GL:KC" << "\t"; // FORMAT

			// determine computed genotype
			pair<int,int> genotype = genotype_likelihoods.get_likeliest_genotype();
			if (ignore_imputed && (nr_unique_kmers == 0)) genotype = {-1,-1};
			if ( (genotype.first != -1) && (genotype.second != -1)) {

				// unique maximum and therefore a likeliest genotype exists
				genotyping_outfile << genotype.first << "/" << genotype.second << ":"; // GT

				// output genotype quality
				genotyping_outfile << genotype_likelihoods.get_genotype_quality(genotype.first, genotype.second) << ":"; // GQ
			} else {
				// genotype could not be determined 
				genotyping_outfile << ".:.:"; // GT:GQ
			}

			// output genotype likelihoods
			vector<long double> likelihoods = genotype_likelihoods.get_all_likelihoods(nr_alleles);
			if (likelihoods.size() < 3) {
				ostringstream oss;
				oss << "Graph::write_genotypes_of: too few likelihoods (" << likelihoods.size() << ") computed for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			ostringstream oss;
			oss << log10(likelihoods[0]);
			for (size_t j = 1; j < likelihoods.size(); ++j) {
				oss << "," << setprecision(4) << log10(likelihoods[j]);
			}
			genotyping_outfile << oss.str(); // GL
			genotyping_outfile << ":" << coverage << endl; // KC
			counter += 1;
		}
	}
}

void Graph::write_phasing(string filename, const vector<GenotypingResult>& genotyping_result, bool write_header, string sample, bool ignore_imputed) {
	if (this->variants_deleted) {
		throw runtime_error("Graph::write_phasing_of: variants have been deleted by delete_variant funtion. Re-build object.");
	}

	if (genotyping_result.size() != this->size()) {
		throw runtime_error("Graph::write_phasing_of: number of variants and number of computed phasings differ.");
	}

	ofstream phasing_outfile;
	if (write_header) {
		phasing_outfile.open(filename);
		if (! phasing_outfile.is_open()) {
			throw runtime_error("Graph::write_phasing_of: phasing output file cannot be opened. Note that the filename must not contain non-existing directories.");
		}

		// write VCF header lines
		phasing_outfile << "##fileformat=VCFv4.2" << endl;
		phasing_outfile << "##fileDate=" << graph_get_date() << endl;
		// TODO output command line
		phasing_outfile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
		phasing_outfile << "##INFO=<ID=UK,Number=1,Type=Integer,Description=\"Total number of unique kmers.\">" << endl;
		phasing_outfile << "##INFO=<ID=AK,Number=R,Type=Integer,Description=\"Number of unique kmers per allele. Will be -1 for alleles not covered by any input haplotype path.\">" << endl;
		phasing_outfile << "##INFO=<ID=MA,Number=1,Type=Integer,Description=\"Number of alleles missing in panel haplotypes.\">" << endl;
		phasing_outfile << "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">" << endl;
		phasing_outfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
		phasing_outfile << "##FORMAT=<ID=KC,Number=1,Type=Float,Description=\"Local kmer coverage.\">" << endl;
		phasing_outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample << endl;
	} else {
		phasing_outfile.open(filename, std::ios_base::app);
		if (! phasing_outfile.is_open()) {
			throw runtime_error("Graph::write_phasing_of: phasing output file cannot be opened. Note that the filename must not contain non-existing directories.");
		}
	}

	size_t counter = 0;
	for (size_t i = 0; i < this->size(); ++i) {
		shared_ptr<Variant> variant = this->variants.at(i);
		unsigned short coverage = genotyping_result.at(i).coverage();
		unsigned short nr_unique_kmers = genotyping_result.at(i).nr_unique_kmers();

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
		variant->separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);

		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			phasing_outfile << v.get_chromosome() << "\t"; // CHROM
			phasing_outfile << (v.get_start_position() + 1) << "\t"; // POS
			phasing_outfile << v.get_id() << "\t"; // ID
			phasing_outfile << v.get_allele_string(0) << "\t"; // REF

			// get alternative alleles
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "Graph::write_phasing_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			vector<string> alt_alleles;
			vector<unsigned char> defined_alleles = {0};
			for (size_t i = 1; i < nr_alleles; ++i) {
				DnaSequence allele = v.get_allele_sequence(i);
				// skip alleles that are undefined
				if (! v.is_undefined_allele(i)) {
					alt_alleles.push_back(allele.to_string());
					defined_alleles.push_back(i);
				}
			}

			string alt_string = "";
			for (unsigned char a = 0; a < alt_alleles.size(); ++a) {
				if (a > 0) alt_string += ',';
				alt_string += alt_alleles[a];
			}

			size_t nr_missing = v.nr_missing_alleles();
			GenotypingResult genotype_likelihoods = singleton_likelihoods.at(j);
			if (nr_missing > 0) genotype_likelihoods = singleton_likelihoods.at(j).get_specific_likelihoods(defined_alleles);

			phasing_outfile << alt_string << "\t"; // ALT
			phasing_outfile << ".\t"; // QUAL
			phasing_outfile << "PASS" << "\t"; // FILTER
			// output allele frequencies of all alleles
			ostringstream info;
			info << "AF=";
			for (unsigned int a = 1; a < defined_alleles.size(); ++a) {
				if (a > 1) info << ",";
				info << setprecision(6) << v.allele_frequency(defined_alleles[a], this->add_reference);				
			}
			info << ";UK=" << nr_unique_kmers; // UK
			info << ";MA=" << nr_missing;

			// if IDs were given in input, write them to output as well
			if (!this->variant_ids[counter].empty()) info << ";ID=" << get_ids(alt_alleles, counter, false);

			phasing_outfile << info.str() << "\t"; // INFO
			phasing_outfile << "GT:KC" << "\t"; // FORMAT

			// determine phasing
			if (ignore_imputed && (nr_unique_kmers == 0)){
				phasing_outfile << "./."; // GT (phased)
			} else {
				pair<unsigned char,unsigned char> haplotype = singleton_likelihoods.at(j).get_haplotype();
				// check if the haplotype allele is undefined
				bool hap1_undefined = v.is_undefined_allele(haplotype.first);
				bool hap2_undefined = v.is_undefined_allele(haplotype.second);
				string hap1 (1, haplotype.first);
				string hap2 (1, haplotype.second);
				if (hap1_undefined) hap1 = ".";
				if (hap2_undefined) hap2 = ".";
				phasing_outfile << (unsigned int) haplotype.first << "|" << (unsigned int) haplotype.second; // GT (phased)
			}
			phasing_outfile << ":" << coverage << endl; // KC
			counter += 1;
		}
	}
}

void Graph::get_left_overhang(size_t index, size_t length, DnaSequence& result) const {
	if (this->variants.at(index) == nullptr) {
		throw runtime_error("Graph::get_left_overhang: variants have been deleted by delete_variant funtion. Re-build object.");
	}

	if (index > 0) {
		if  (this->variants.at(index-1) == nullptr) {
			throw runtime_error("Graph::get_left_overhang: variants have been deleted by delete_variant funtion. Re-build object.");
		}
	}

	size_t cur_start = this->variants.at(index)->get_start_position();
	size_t prev_end = 0;
	if (index > 0) prev_end = this->variants.at(index-1)->get_end_position();
	size_t overhang_start = cur_start - length;
	if (overhang_start < prev_end) overhang_start = prev_end;
	size_t overhang_end = cur_start;
	this->fasta_reader.get_subsequence(this->chromosome, overhang_start, overhang_end, result);
}

void Graph::get_right_overhang(size_t index, size_t length, DnaSequence& result) const {
	if (this->variants.at(index) == nullptr) {
		throw runtime_error("Graph::get_right_overhang: variants have been deleted by delete_variant funtion. Re-build object.");
	}

	if (index < (this->variants.size() - 1)) {
		if (this->variants.at(index+1) == nullptr){
			throw runtime_error("Graph::get_right_overhang: variants have been deleted by delete_variant funtion. Re-build object.");
		}
	}

	size_t cur_end = this->variants.at(index)->get_end_position();
	size_t next_start = this->fasta_reader.get_size_of(this->chromosome);
	if (index < (this->variants.size() - 1)) next_start = this->variants.at(index+1)->get_start_position();
	size_t overhang_end = cur_end + length;
	if (overhang_end > next_start) overhang_end = next_start;
	size_t overhang_start = cur_end;
	this->fasta_reader.get_subsequence(this->chromosome, overhang_start, overhang_end, result);
}

void Graph::delete_variant(size_t index) {
	if (index < this->size()) {
		if (this->variants.at(index) != nullptr) {
			assert(this->variants[index].unique());
			this->variants[index].reset();
			this->variants[index] = nullptr;
			this->variants_deleted = true;
		}
	} else {
		throw runtime_error("Graph::get_variant: index out of bounds.");
	}
}

bool Graph::variants_were_deleted() const {
	return this->variants_deleted;
}