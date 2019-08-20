#include <sstream>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <math.h>
#include "variantreader.hpp"

using namespace std;

void parse_line(vector<DnaSequence>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(DnaSequence(token));
	}
}

void parse_line(vector<string>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

VariantReader::VariantReader(string filename, string reference_filename, size_t kmer_size, string sample)
	:fasta_reader(reference_filename),
	 kmer_size(kmer_size),
	 nr_variants(0),
	 sample(sample),
	 genotyping_outfile_open(false),
	 phasing_outfile_open(false)
{
	if (filename.substr(filename.size()-3,3).compare(".gz") == 0) {
		throw runtime_error("VariantReader::VariantReader: Uncompressed VCF-file is required.");
	}
	ifstream file(filename);
	string line;
	string previous_chrom("");
	size_t previous_end_pos = 0;
	map<unsigned int, string> fields = { {0, "#CHROM"}, {1, "POS"}, {2, "ID"}, {3, "REF"}, {4, "ALT"}, {5, "QUAL"}, {6, "FILTER"}, {7, "INFO"}, {8, "FORMAT"} };
	vector<Variant> variant_cluster;
	// read VCF-file line by line
	while (getline(file, line)) {
		if (line.size() == 0) continue;
		vector<string> tokens;
		parse_line(tokens, line, '\t');
		// header lines?
		if (tokens[0].substr(0,2) == "##") continue;
		if (tokens[0].at(0) == '#') {
			// check number of samples/paths given
			if (tokens.size() < 9) {
				throw runtime_error("VariantReader::VariantReader: not a proper VCF-file.");
			}
			if (tokens.size() < 10) {
				throw runtime_error("VariantReader::VariantReader: no haplotype paths given.");
			}
			// validate header line
			for (unsigned int i = 0; i < 9; ++i) {
				if (tokens[i] != fields[i]) {
					throw runtime_error("VariantReader::VariantReader: VCF header line is malformed.");
				}
			}
			this->nr_paths = tokens.size() - 9;
			continue;
		}
		if (tokens.size() < 10) {
			throw runtime_error("VariantReader::VariantReader: malformed VCF-file, or no haplotype paths given in VCF.");
		}
		// get chromosome
		string current_chrom = tokens[0];
		// get position
		size_t current_start_pos;
		stringstream sstream(tokens[1]);
		sstream >> current_start_pos;
		// VCF positions are 1-based
		current_start_pos -= 1;
		// if distance to next variant is larger than kmer_size, start a new cluster
		if ( (previous_chrom != current_chrom) || (current_start_pos - previous_end_pos) >= (kmer_size-1) ) {
			// merge all variants currently in cluster and store them
			add_variant_cluster(previous_chrom, &variant_cluster);
			variant_cluster.clear();
		}
		// get REF allele
		DnaSequence ref(tokens[3]);
		DnaSequence observed_allele;
		this->fasta_reader.get_subsequence(current_chrom, current_start_pos, current_start_pos + ref.size(), observed_allele);
		if (ref != observed_allele) {
			throw runtime_error("VariantReader::VariantReader: reference allele given in VCF does not match allele in reference fasta file at that position.");
		}
		size_t current_end_pos = current_start_pos + ref.size();
		// get ALT alleles
		vector<DnaSequence> alleles = {ref};
		parse_line(alleles, tokens[4], ',');
		// construct paths
		vector<unsigned char> paths;
		for (size_t i = 9; i < tokens.size(); ++i) {
			// make sure all genotypes are phased
			if (tokens[i].find('/') != string::npos) {
				throw runtime_error("VariantReader::VariantReader: Found unphased genotype.");
			}
			vector<string> p;
			parse_line(p, tokens[i], '|');
			for (string& s : p){
				paths.push_back( (unsigned char) atoi(s.c_str()));
			}
		}
		// determine left and right flanks
		DnaSequence left_flank;
		this->fasta_reader.get_subsequence(current_chrom, current_start_pos - kmer_size + 1, current_start_pos, left_flank);
		DnaSequence right_flank;
		this->fasta_reader.get_subsequence(current_chrom, current_end_pos, current_end_pos + kmer_size - 1, right_flank);
		// add Variant to variant_cluster
		Variant variant (left_flank, right_flank, current_chrom, current_start_pos, current_end_pos, alleles, paths);
		variant_cluster.push_back(variant);
		previous_chrom = current_chrom;
		previous_end_pos = current_end_pos;
	}
	// add last cluster to list
	add_variant_cluster(previous_chrom, &variant_cluster);
	cerr << "Identified " << this->nr_variants << " variants in total from VCF-file." << endl;
}

size_t VariantReader::get_kmer_size() const {
	return this->kmer_size;
}

void VariantReader::write_path_segments(std::string filename) const {
	ofstream outfile;
	outfile.open(filename);
	for (auto element : this->variants_per_chromosome) {
		assert(element.second.size() > 0);
		size_t prev_end = 0;
		for (auto variant : element.second) {
			// generate reference unitig and write to file
			size_t start_pos = variant.get_start_position();
			outfile << ">" << element.first << "_reference_" << start_pos << endl;
			string ref_segment;
			this->fasta_reader.get_subsequence(element.first, prev_end, start_pos, ref_segment);
			outfile << ref_segment << endl;
			for (size_t allele = 0; allele < variant.nr_of_alleles(); ++allele) {
				// sequence name
				outfile << ">" << element.first << "_" << start_pos << "_" << allele << endl;
				outfile << variant.get_allele_string(allele) << endl;
			}
			prev_end = variant.get_end_position();
		}
		// output reference sequence after last position on chromosome
		outfile << ">" << element.first << "_reference_end" << endl;
		size_t chr_len = this->fasta_reader.get_size_of(element.first);
		string ref_segment;
		this->fasta_reader.get_subsequence(element.first, prev_end, chr_len, ref_segment);
		outfile << ref_segment << endl;
	}
	outfile.close();
}

void VariantReader::get_chromosomes(vector<string>* result) const {
	for (auto const& element : this->variants_per_chromosome) {
		result->push_back(element.first);
	}
}

size_t VariantReader::size_of(string chromosome) const {
	if (this->variants_per_chromosome.find(chromosome) != this->variants_per_chromosome.end()) {
		return this->variants_per_chromosome.at(chromosome).size();
	} else {
		return 0;
	}
}

const Variant& VariantReader::get_variant(string chromosome, size_t index) const {
	if (index < size_of(chromosome)) {
		return this->variants_per_chromosome.at(chromosome)[index];
	} else {
		throw runtime_error("VariantReader::get_variant: index out of bounds.");
	}
}

void VariantReader::add_variant_cluster(string& chromosome, vector<Variant>* cluster) {
	if (!cluster->empty()) {
		// merge all variants in cluster
		Variant combined = cluster->at(0);
		for (size_t v = 1; v < cluster->size(); ++v) {
			combined.combine_variants(cluster->at(v));
		}
		combined.add_flanking_sequence();
		this->variants_per_chromosome[chromosome].push_back(combined);
		this->nr_variants += 1;
	}
}

const vector<Variant>& VariantReader::get_variants_on_chromosome(string chromosome) const {
	if (this->variants_per_chromosome.find(chromosome) != this->variants_per_chromosome.end()) {
		return this->variants_per_chromosome.at(chromosome);
	} else {
		throw runtime_error("VariantReader::get_variants_on_chromosome: chromosome " + chromosome + " not present in VCF.");
	}
}

string get_date() {
	time_t t = time(0);
	tm* now = localtime(&t);
	ostringstream oss;
	oss << (now->tm_year+1900) << setfill('0') << setw(2) << (now->tm_mon+1) << setw(2) << now->tm_mday;
	return string(oss.str());
}

void VariantReader::open_genotyping_outfile(string filename) {
	this->genotyping_outfile.open(filename);
	if (! this->genotyping_outfile.is_open()) {
		throw runtime_error("VariantReader::open_genotyping_outfile: genotyping output file cannot be opened.");
	}

	this->genotyping_outfile_open = true;
	
	// write VCF header lines
	this->genotyping_outfile << "##fileformat=VCFv4.2" << endl;
	this->genotyping_outfile << "##fileDate=" << get_date() << endl;
	// TODO output command line
	this->genotyping_outfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	this->genotyping_outfile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality: phred scaled probability that the genotype is wrong.\">" << endl;
	this->genotyping_outfile << "##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Comma-separated log10-scaled genotype likelihoods for absent, heterozygous, homozygous.\">" << endl;
	this->genotyping_outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << this->sample << endl; 
}

void VariantReader::open_phasing_outfile(string filename) { 
	this->phasing_outfile.open(filename);
	if (! this->phasing_outfile.is_open()) {
		throw runtime_error("VariantReader::open_phasing_outfile: phasing output file cannot be opened.");
	}

	this->phasing_outfile_open = true;

	// write VCF header lines
	this->phasing_outfile << "##fileformat=VCFv4.2" << endl;
	this->phasing_outfile << "##fileDate=" << get_date() << endl;
	// TODO output command line
	this->phasing_outfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	this->phasing_outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << this->sample << endl;
}

void VariantReader::write_genotypes_of(string chromosome, const vector<GenotypingResult>& genotyping_result) {
	// outfile needs to be open
	if (!this->genotyping_outfile_open) {
		throw runtime_error("VariantReader::write_genotypes_of: output file needs to be opened before writing.");
	}

	if (this->variants_per_chromosome.find(chromosome) == this->variants_per_chromosome.end()) {
		cerr << "VariantReader::write_genotypes_of: no variants for given chromosome were written." << endl;
		return;
	}

	if (genotyping_result.size() != size_of(chromosome)) {
		throw runtime_error("VariantReader::write_genotypes_of: number of variants and number of computed genotypes differ.");
	}
	
	for (size_t i = 0; i < size_of(chromosome); ++i) {
		Variant variant = this->variants_per_chromosome.at(chromosome)[i];

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
		variant.separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);

		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			this->genotyping_outfile << v.get_chromosome() << "\t"; // CHROM
			this->genotyping_outfile << (v.get_start_position() + 1) << "\t"; // POS
			this->genotyping_outfile << ".\t"; // ID
			this->genotyping_outfile << v.get_allele_string(0) << "\t"; // REF

			// get alternative alleles
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "VariantReader::write_genotypes_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			string alt_alleles = v.get_allele_string(1);
			for (size_t i = 2; i < nr_alleles; ++i) {
				alt_alleles += "," + v.get_allele_string(i);
			}
			this->genotyping_outfile << alt_alleles << "\t"; // ALT
			this->genotyping_outfile << ".\t"; // QUAL
			this->genotyping_outfile << "PASS" << "\t"; // FILTER
			this->genotyping_outfile << "." << "\t"; // INFO
			this->genotyping_outfile << "GT:GQ:GL" << "\t"; // FORMAT

			// determine computed genotype
			pair<int,int> genotype = singleton_likelihoods.at(j).get_likeliest_genotype();
			if ( (genotype.first != -1) && (genotype.second != -1)) {

				// unique maximum and therefore a likeliest genotype exists
				this->genotyping_outfile << genotype.first << "/" << genotype.second << ":"; // GT

				// output genotype quality
				this->genotyping_outfile << singleton_likelihoods.at(j).get_genotype_quality(genotype.first, genotype.second) << ":"; // GQ
			} else {
				// genotype could not be determined 
				this->genotyping_outfile << ".:.:"; // GT:GQ
			}

			// output genotype likelihoods
			vector<long double> likelihoods = singleton_likelihoods.at(j).get_all_likelihoods(nr_alleles);
			if (likelihoods.size() < 3) {
				ostringstream oss;
				oss << "VariantReader::write_genotypes_of: too few likelihoods (" << likelihoods.size() << ") computed for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			ostringstream oss;
			oss << log10(likelihoods[0]);
			for (size_t j = 1; j < likelihoods.size(); ++j) {
				oss << "," << setprecision(4) << log10(likelihoods[j]);
			}
			this->genotyping_outfile << oss.str() << endl; // GL
		}
	}
}

void VariantReader::write_phasing_of(string chromosome, const vector<GenotypingResult>& genotyping_result) {
	// outfile needs to be open
	if (! this->phasing_outfile_open) {
		throw runtime_error("VariantReader::write_phasing_of: output file needs to be opened before writing.");
	}
	if (genotyping_result.size() != size_of(chromosome)) {
		throw runtime_error("VariantReader::write_phasing_of: number of variants and number of computed phasings differ.");
	}

	for (size_t i = 0; i < size_of(chromosome); ++i) {
		Variant variant = this->variants_per_chromosome.at(chromosome)[i];

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
		variant.separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);

		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			this->phasing_outfile << v.get_chromosome() << "\t"; // CHROM
			this->phasing_outfile << (v.get_start_position() + 1) << "\t"; // POS
			this->phasing_outfile << ".\t"; // ID
			this->phasing_outfile << v.get_allele_string(0) << "\t"; // REF

			// get alternative alleles
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "VariantReader::write_phasing_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			string alt_alleles = v.get_allele_string(1);
			for (size_t i = 2; i < nr_alleles; ++i) {
				alt_alleles += "," + v.get_allele_string(i);
			}

			this->phasing_outfile << alt_alleles << "\t"; // ALT
			this->phasing_outfile << ".\t"; // QUAL
			this->phasing_outfile << "PASS" << "\t"; // FILTER
			this->phasing_outfile << "." << "\t"; // INFO
			this->phasing_outfile << "GT" << "\t"; // FORMAT

			// determine phasing
			pair<unsigned char,unsigned char> haplotype = singleton_likelihoods.at(j).get_haplotype();
			this->phasing_outfile << (unsigned int) haplotype.first << "|" << (unsigned int) haplotype.second << endl; // GT (phased)
		}
	}
}

void VariantReader::close_genotyping_outfile() {
	this->genotyping_outfile.close();
	this->genotyping_outfile_open = false;
}

void VariantReader::close_phasing_outfile() {
	this->phasing_outfile.close();
	this->phasing_outfile_open = false;
}

size_t VariantReader::nr_of_genomic_kmers() const {
	return this->fasta_reader.get_total_kmers(this->kmer_size);
}
