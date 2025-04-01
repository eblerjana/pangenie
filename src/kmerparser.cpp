#include<iostream>
#include<cassert>
#include <sstream>
#include "kmerparser.hpp"

using namespace std;

void parse(vector<string>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

void parse_kmer_line(string line, string& chrom, size_t& start, vector<string>& kmers, vector<string>& flanking_kmers, bool& is_header) {
	vector<string> tokens;
	parse(tokens, line, '\t');
	assert (tokens.size() == 5);
	if (tokens[0][0] == '#') {
		is_header = true;
		return;
	}
	chrom = tokens[0];
	start = atoi(tokens[1].c_str());
	if (tokens[3] != "nan") parse(kmers, tokens[3], ',');
	if (tokens[4] != "nan") parse(flanking_kmers, tokens[4], ',');
}

unsigned short compute_local_coverage(vector<string>& kmers, shared_ptr<KmerCounter> read_counts, size_t kmer_coverage) {
	size_t total_coverage = 0;
	size_t total_kmers = 0;
	size_t min_cov = kmer_coverage / 4;
	size_t max_cov = kmer_coverage * 4;

	for (auto& kmer : kmers) {
		size_t read_count = read_counts->getKmerAbundance(kmer);
		// ignore too extreme counts
		if ( (read_count < min_cov) || (read_count > max_cov) ) continue;
		total_coverage += read_count;
		total_kmers += 1;		
	}

	if ((total_kmers > 0) && (total_coverage > 0)){
		return total_coverage / total_kmers;
	} else {
		return kmer_coverage;
	}		
}