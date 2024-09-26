#include "catch.hpp"
#include "utils.hpp"
#include "../src/commands.hpp"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

TEST_CASE("Commands run_genotype_command") {
	string precomputed_prefix = "TODO"
	string readfile = "TODO"
	string outname = "test"
	string sample_name = "sample"
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.01L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	size_t sampling_size = 0;
	uint64_t hash_size = 3000000000;
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
	REQUIRE(computed_lines.size() == 4);
	vector<string> expected_line1;
	vector<string> expected_line2;
	vector<string> expected_line3;
	vector<string> expected_line4;

	REQUIRE(computed_lines[0] == expected_lines1);
	REQUIRE(computed_lines[1] == expected_lines2);
	REQUIRE(computed_lines[2] == expected_lines3);
	REQUIRE(computed_lines[3] == expected_lines4);

}