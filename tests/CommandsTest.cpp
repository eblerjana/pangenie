#include "catch.hpp"
#include "utils.hpp"
#include "../src/commands.hpp"
#include "../src/probabilitytable.hpp"
#include "../src/hmm.hpp"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cereal/archives/binary.hpp>

using namespace std;

TEST_CASE("Commands run_genotype_command") {
	string precomputed_prefix = "../tests/data/index";
	string readfile = "../tests/data/region-reads.fa";
	string outname = "../tests/data/testfull";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 0;
	uint64_t hash_size = 100000;
	size_t panel_size = 0;
	double recombrate = 1.26;
	bool output_panel = false;

	run_genotype_command(precomputed_prefix, readfile, outname, sample_name, nr_jellyfish_threads, nr_core_threads, only_genotyping, only_phasing, effective_N, regularization, count_only_graph, ignore_imputed, sampling_size, hash_size, panel_size, recombrate, output_panel);

	// check if output file exists
	{
		ifstream file(outname + "_genotyping.vcf");
		REQUIRE(file.good());
	}

	// parse output file
	string line;
	vector<vector<string>> computed_lines;
	parse_vcf_lines(outname + "_genotyping.vcf");

	// check if output looks as expected
	REQUIRE(computed_lines.size() == 2);


	REQUIRE(computed_lines[0] == expected_lines1);
	REQUIRE(computed_lines[1] == expected_lines2);
	REQUIRE(computed_lines[2] == expected_lines3);
	REQUIRE(computed_lines[3] == expected_lines4);


	// Determine results to expect
	size_t kmer_abundance_peak = 18;
	ProbabilityTable probs = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);
	UniqueKmersMap uk;

	// Reconstruct expected UniqueKmers objects from file
	ifstream is("TEST_UniqueKmersList.cereal", std::ios::binary);
	cereal::BinaryInputArchive archive_is( is );
	archive_is(uk);

	HMM hmm(&uk.unique_kmers["chr1"], &probs, only_genotyping, only_phasing, recombrate, false, effective_N);
	vector<GenotypingResult> genotypes = hmm.get_genotyping_result();

	REQUIRE(genotypes.size() == 2);
	genotypes[0].normalize();
	genotypes[1].normalize();
	vector<unsigned char> defined1 = {0,1};
	vector<unsigned char> defined2 = {0,1,2};
	vector<vector<unsigned char>> defined = {{0,1}, {0,1,2}};
	cout << "COMMMMMMM" << endl;

	vector<string> expected_likelihoods = {};

	for(size_t i = 0; i < 2; ++i) {
		vector<long double> likelihoods = genotypes[i].get_specific_likelihoods(defined[i]).get_all_likelihoods(defined[i].size());
		for (auto g : likelihoods) {
			ostringstream oss;
			oss << setprecision(4) << log10(g) << endl;
			expected_likelihoods.push_back(oss.str());
		}
	}





}