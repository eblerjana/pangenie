#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <regex>
#include "graphbuilder.hpp"

using namespace std;

void builder_parse_line(vector<DnaSequence>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(DnaSequence(token));
	}
}

void builder_parse_line(vector<string>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

bool builder_matches_pattern(string& sequence, regex& r) {
	size_t total_length = sequence.size();
	string control = "";
	auto begin = sequence.begin();
	size_t length = 1000;
	size_t letters_seen = 0;
	while (begin + letters_seen != sequence.end()) {
		size_t next_interval = min(length, total_length - letters_seen);
		if (!regex_match(begin + letters_seen, begin + letters_seen + next_interval, r)) {
			return false;
		}
		control += string(begin + letters_seen, begin + letters_seen + next_interval);
		letters_seen += next_interval;
	}
	assert(control == sequence);
	return true;
}

void builder_parse_info_fields(vector<string>& result, string line) {
	vector<string> fields;
	builder_parse_line(fields, line, ';');
	for (auto s : fields) {
		if (s.rfind("ID=", 0) == 0) {
			// ID field present
			builder_parse_line(result, s.substr (3), ',');
		}
	}
}

GraphBuilder::GraphBuilder(string filename, string reference_filename, map<string, shared_ptr<Graph>>& result, string segments_file,  size_t kmer_size, bool add_reference)
	: kmer_size(kmer_size),
	  nr_variants(0)
{
	cerr << "Read reference genome ..." << endl;
	// read the reference sequence
	FastaReader fasta_reader(reference_filename);
	// read variants from input VCF and merge such closer than kmer size into bubbles
	cerr << "Read input VCF ..." << endl;
	construct_graph(filename, &fasta_reader, result, add_reference);
	cerr << "Write path segments to file ..." << endl;
	// write a FASTA containing all graph sequences (needed for kmer counting)
	write_path_segments(segments_file, &fasta_reader, result);
}

void GraphBuilder::construct_graph(std::string filename, FastaReader* fasta_reader, std::map<string, shared_ptr<Graph>>& result, bool add_reference)
{
	// stores chromosome names and their sizes (= nr of variant bubbles)
	vector<pair<size_t,string>> chromosome_sizes;

	if (filename.substr(filename.size()-3,3).compare(".gz") == 0) {
		throw runtime_error("GraphBuilder::GraphBuilder: Uncompressed VCF-file is required.");
	}
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("GraphBuilder::GraphBuilder: input VCF file cannot be opened.");
	}
	string line;
	string previous_chrom("");
	size_t previous_end_pos = 0;
	map<unsigned int, string> fields = { {0, "#CHROM"}, {1, "POS"}, {2, "ID"}, {3, "REF"}, {4, "ALT"}, {5, "QUAL"}, {6, "FILTER"}, {7, "INFO"}, {8, "FORMAT"} };

	// variants to be merged into a bubble (since they are less than the kmersize apart)
	vector<shared_ptr<Variant>> variant_cluster;
	// IDs of individual variants alleles
	vector<vector<string>> variant_cluster_ids;

	// Graph object of the current chromosome
	shared_ptr<Graph> current_graph = nullptr;

	// read VCF-file line by line
	while (getline(file, line)) {
		if (line.size() == 0) continue;
		vector<string> tokens;
		builder_parse_line(tokens, line, '\t');
		// header lines?
		if (tokens[0].substr(0,2) == "##") continue;
		if (tokens[0].at(0) == '#') {
			// check number of samples/paths given
			if (tokens.size() < 9) {
				throw runtime_error("GraphBuilder::GraphBuilder: not a proper VCF-file.");
			}
			if (tokens.size() < 10) {
				throw runtime_error("GraphBuilder::GraphBuilder: no haplotype paths given.");
			}
			// validate header line
			for (unsigned int i = 0; i < 9; ++i) {
				if (tokens[i] != fields[i]) {
					throw runtime_error("GraphBuilder::GraphBuilder: VCF header line is malformed.");
				}
			}
			this->nr_paths = (tokens.size() - 9)*2;
			// add one for reference path
			if (add_reference) this->nr_paths += 1;
			continue;
		}
		if (tokens.size() < 10) {
			throw runtime_error("GraphBuilder::GraphBuilder: malformed VCF-file, or no haplotype paths given in VCF.");
		}
		// get chromosome
		string current_chrom = tokens[0];
		// get position
		size_t current_start_pos;
		stringstream sstream(tokens[1]);
		sstream >> current_start_pos;
		// VCF positions are 1-based
		current_start_pos -= 1;
		// if variant is contained in previous one, skip it
		if ((previous_chrom == current_chrom) && (current_start_pos < previous_end_pos)) {
			stringstream err_msg;
			err_msg << "GraphBuilder: variant at " << current_chrom << ":" << current_start_pos << " overlaps previous one. VCF does not represent a pangenome graph."  << endl;
			throw runtime_error(err_msg.str());
		}

		// get REF allele
		DnaSequence ref(tokens[3]);
		DnaSequence observed_allele;

		if (previous_chrom == current_chrom) {
			current_graph->get_fasta_reader().get_subsequence(current_chrom, current_start_pos, current_start_pos + ref.size(), observed_allele);
		} else {
			fasta_reader->get_subsequence(current_chrom, current_start_pos, current_start_pos + ref.size(), observed_allele);
		}

		if (ref != observed_allele) {
			throw runtime_error("GraphBuilder::GraphBuilder: reference allele given in VCF does not match allele in reference fasta file at that position.");
		}
		size_t current_end_pos = current_start_pos + ref.size();
		// get ALT alleles
		vector<DnaSequence> alleles = {ref};
		// make sure alt alleles are given explicitly
		regex r("^[CAGTcagt,]+$");
		if (!builder_matches_pattern(tokens[4], r)) {
			// skip this position
			cerr << "GraphBuilder: skip variant at " << current_chrom << ":" << current_start_pos << " since alleles contain undefined nucleotides: " << tokens[4] << endl;
			continue;
		}
		builder_parse_line(alleles, tokens[4], ',');

		// currently, number of alleles is limited to 256
		if (alleles.size() > 255) {
			throw runtime_error("GraphBuilder: number of alternative alleles is limited to 254 in current implementation. Make sure the VCF contains only alternative alleles covered by at least one of the haplotypes.");
		}

		// determine size of current chromosome. If Graph object was created already, FastaReader no longer contains
		// the corresponding sequence and we have to get the info from the Graph
		size_t size_of_chromosome = 0;
		if (previous_chrom == current_chrom) {
			size_of_chromosome = current_graph->get_fasta_reader().get_size_of(current_chrom);
		} else {
			size_of_chromosome = fasta_reader->get_size_of(current_chrom);
		}

		// TODO: handle cases where variant is less than kmersize from start or end of the chromosome
		if ( (current_start_pos < (kmer_size*2) ) || ( (current_end_pos + (kmer_size*2)) > size_of_chromosome) ) {
			cerr << "GraphBuilder: skip variant at " << current_chrom << ":" << current_start_pos << " since variant is less than 2 * kmer size from start or end of chromosome. " << endl;

			continue;
		}

		// if distance to next variant is larger than kmer_size or the chromosome changed, start a new cluster
		if ( (previous_chrom != current_chrom) || (current_start_pos - previous_end_pos) >= (this->kmer_size-1) ) {
			// merge all variants currently in cluster and store them
			if (current_graph != nullptr) current_graph->add_variant_cluster(&variant_cluster, variant_cluster_ids, true);
			variant_cluster.clear();
			variant_cluster_ids.clear();

			if (previous_chrom != current_chrom) {
				// chromosome changed, construct a Graph object for new chromosome
				if (current_graph != nullptr) {
					result[previous_chrom] = current_graph;
				}
				current_graph = shared_ptr<Graph>(new Graph(fasta_reader->extract_name(current_chrom), current_chrom, this->kmer_size, add_reference));
			}
		}

		// store mapping of alleles to variant ids
		vector<string> var_ids;
		builder_parse_info_fields(var_ids, tokens[7]);

		// make sure that there are at most 255 paths (including reference path in case it is requested)
		if (this->nr_paths > 255) {
			throw runtime_error("GraphBuilder: number of paths is limited to 254 in current implementation.");
		}

		// construct paths
		vector<unsigned char> paths = {};
		if (add_reference) paths.push_back((unsigned char) 0);
		unsigned char undefined_index = alleles.size();
		string undefined_allele = "N";

		for (size_t i = 9; i < tokens.size(); ++i) {
			// make sure all genotypes are phased
			if (tokens[i].find('/') != string::npos) {
				throw runtime_error("GraphBuilder::GraphBuilder: Found unphased genotype.");
			}
			vector<string> p ;
			builder_parse_line(p, tokens[i], '|');
			if (p.size() != 2) {
				throw runtime_error("GraphBuilder::GraphBuilder: Found invalid genotype. Genotypes must be diploid (.|. if missing).");
			}
			for (string& s : p){
				// handle unknown genotypes '.'
				if (s == ".") {
					// add "N" allele to the list of alleles
					builder_parse_line(alleles, undefined_allele, ',');
					paths.push_back(undefined_index);
					assert(undefined_index < 255);
					undefined_index += 1;
				} else {
					unsigned int p_index = atoi(s.c_str());
					if (p_index >= alleles.size()) {
						throw runtime_error("GraphBuilder::GraphBuilder: invalid genotype in VCF.");
					}
					assert(p_index < 255);
					paths.push_back( (unsigned char) p_index);
				}
			}
		}

		assert(current_graph != nullptr);

		// determine left and right flanks
		DnaSequence left_flank;
		current_graph->get_fasta_reader().get_subsequence(current_chrom, current_start_pos - kmer_size + 1, current_start_pos, left_flank);
		DnaSequence right_flank;
		current_graph->get_fasta_reader().get_subsequence(current_chrom, current_end_pos, current_end_pos + kmer_size - 1, right_flank);
		// add Variant to variant_cluster
		shared_ptr<Variant> variant = shared_ptr<Variant>(new Variant(left_flank, right_flank, current_chrom, current_start_pos, current_end_pos, alleles, paths));
		variant_cluster.push_back(variant);
		variant_cluster_ids.push_back(var_ids);
		previous_chrom = current_chrom;
		previous_end_pos = current_end_pos;

	}

	// add last cluster to list and store the Graph object
	if (current_graph != nullptr) {
		current_graph->add_variant_cluster(&variant_cluster, variant_cluster_ids, true);
		result[previous_chrom] = current_graph;
	}

	// determine total number of variant clusters read and store chromosomes in order of their size
	for (auto it = result.begin(); it != result.end(); ++it) {
		chromosome_sizes.push_back(make_pair(it->second->size(), it->first));
		this->nr_variants += it->second->size();
	}

	sort(chromosome_sizes.rbegin(), chromosome_sizes.rend());
	for (auto const& element : chromosome_sizes) {
		this->chromosomes.push_back(element.second);
	}

	cerr << "Identified " << this->nr_variants << " variants in total from VCF-file." << endl;
}

size_t GraphBuilder::get_kmer_size() const {
	return this->kmer_size;
}

void GraphBuilder::get_chromosomes(vector<string>* result) const {
	for (auto c : this->chromosomes) result->push_back(c);
}

size_t GraphBuilder::nr_of_paths() const {
	return this->nr_paths;
}

void GraphBuilder::write_path_segments(string filename, FastaReader* fasta_reader, map<string, shared_ptr<Graph>>& result) const {
	ofstream outfile;
	outfile.open(filename);
	if (!outfile.good()) {
		stringstream ss;
		ss << "GraphBuilder::write_path_segments: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
		throw runtime_error(ss.str());
	}

	// collect all chromosomes
	vector<string> vcf_chromosomes;
	this->get_chromosomes(&vcf_chromosomes); 

	// make sure to capture all chromosomes in the reference (including such for which no variants are given)
	vector<string> all_chromosome_names = vcf_chromosomes;
	fasta_reader->get_sequence_names(all_chromosome_names);

	for (auto element : all_chromosome_names) {
		size_t prev_end = 0;
		// check if chromosome was present in VCF and write allele sequences in this case
		if (find(vcf_chromosomes.begin(), vcf_chromosomes.end(), element) != vcf_chromosomes.end()) {
			if (result.at(element)->variants_were_deleted()) {
				throw runtime_error("GraphBuilder::write_path_segments: variants have been deleted by delete_variant funtion. Re-build object.");
			}

			const FastaReader& graph_fasta_reader = result.at(element)->get_fasta_reader();
			for (size_t i = 0; i < result.at(element)->size(); ++i) {
				const Variant& variant = result.at(element)->get_variant(i);
				// generate reference unitig and write to file
				size_t start_pos = variant.get_start_position();
				outfile << ">" << element << "_reference_" << start_pos << endl;
				string ref_segment;
				graph_fasta_reader.get_subsequence(element, prev_end, start_pos, ref_segment);
				outfile << ref_segment << endl;
				for (size_t allele = 0; allele < variant.nr_of_alleles(); ++allele) {
					// sequence name
					outfile << ">" << element << "_" << start_pos << "_" << allele << endl;
					outfile << variant.get_allele_string(allele) << endl;
				}
				prev_end = variant.get_end_position();
			}

			// output reference sequence after last position on chromosome
			outfile << ">" << element << "_reference_end" << endl;
			size_t chr_len = graph_fasta_reader.get_size_of(element);
			string ref_segment;
			graph_fasta_reader.get_subsequence(element, prev_end, chr_len, ref_segment);
			outfile << ref_segment << endl;

		} else {
			// output chromosomes not present in VCF
			outfile << ">" << element << "_reference_end" << endl;
			size_t chr_len = fasta_reader->get_size_of(element);
			string ref_segment;
			assert(fasta_reader->contains_name(element));
			fasta_reader->get_subsequence(element, prev_end, chr_len, ref_segment);
			outfile << ref_segment << endl;
		}
	}
	outfile.close();
}