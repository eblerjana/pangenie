#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include "../src/hmm.hpp"
#include "../src/probabilitycomputer.hpp"
#include "../src/genotypingresult.hpp"
#include "utils.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

unsigned char string_to_u_char (string s) {
	istringstream reader (s);
	unsigned char result;
	reader >> result;
	return result;
}

long double string_to_l_double (string s) {
	istringstream reader (s);
	long double result;
	reader >> result;
	return result;
}

TEST_CASE("Posteriors compare_posteriors", "[Posteriors compare_posteriors]") {
	vector<unsigned char> path_to_allele = {0, 1};
	UniqueKmers u1(2000, path_to_allele);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);
	u1.set_coverage(5);

	UniqueKmers u2(3000, path_to_allele);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(5, a2);
	u2.set_coverage(5);

	ProbabilityTable probs(5, 10, 30, 0.0L);
	probs.modify_probability(5, 10, CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(5, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(5, 5, CopyNumber(0.9,0.3,0.1));
	vector<UniqueKmers*> unique_kmers = {&u1,&u2};

	map<string, bool> posterior_pos = { {"chr1_2001", true}, {"chr1_3001", true}};

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25, nullptr, true, &posterior_pos, "test1", "chr1");

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
	}
	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );

	// make sure genotype likelihoods are consistent with the posteriors
	ifstream post_file("test1_chr1_posteriors.tsv");
	REQUIRE(post_file.good());

	string prev_var = "";
	GenotypingResult genotypes;
	bool cleared = false;
	computed_likelihoods.clear();

	string line;
	size_t index = 0;
	while(getline(post_file, line)){
		vector<string> fields;
		string token;
		istringstream iss (line);
		while(getline(iss, token, '\t')) {
			fields.push_back(token);
		}

		REQUIRE(fields.size() == 6);

		// if variant has changed, compare likelihoods
		if ( (fields[0] != prev_var) && (index != 0)) {
			// compare computed likelihoods to expected ones
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 0));
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 1));
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 1, (unsigned char) 1));
			cleared = true;
			genotypes = GenotypingResult();
		}

		// determine which genotype the current state refers to
		unsigned char allele1 = (unsigned char) atoi(fields[4].c_str());
		unsigned char allele2 = (unsigned char) atoi(fields[5].c_str());
		genotypes.add_to_likelihood(allele1, allele2, string_to_l_double(fields[3]));
		cleared = false;
		index += 1;
		prev_var = fields[0];		
	
	}

	if (!cleared) {
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 0));
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 1));
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 1, (unsigned char) 1));
	}

	// variants are enumerated in reverse order, that is why likelihoods also appear in reverse order (v2, v1)
	vector<double> expected_likelihoods_rev = { 0.9678020017, 0.031003181, 0.0011948172, 0.0509465435, 0.9483202731, 0.0007331832};

	REQUIRE(computed_likelihoods.size() == expected_likelihoods_rev.size());
	REQUIRE(compare_vectors(expected_likelihoods_rev, computed_likelihoods));
}

TEST_CASE("Posteriors compare_posteriors_subset", "[Posteriors compare_posteriors_subset]") {
	vector<unsigned char> path_to_allele = {0, 1};
	UniqueKmers u1(2000, path_to_allele);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);
	u1.set_coverage(5);

	UniqueKmers u2(3000, path_to_allele);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(5, a2);
	u2.set_coverage(5);

	ProbabilityTable probs(5, 10, 30, 0.0L);
	probs.modify_probability(5, 10, CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(5, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(5, 5, CopyNumber(0.9,0.3,0.1));
	vector<UniqueKmers*> unique_kmers = {&u1,&u2};

	map<string, bool> posterior_pos = {{"chr1_3001", true}};

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25, nullptr, true, &posterior_pos, "test1", "chr1");

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
	}
	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );

	// make sure genotype likelihoods are consistent with the posteriors
	ifstream post_file("test1_chr1_posteriors.tsv");
	REQUIRE(post_file.good());

	string prev_var = "";
	GenotypingResult genotypes;
	bool cleared = false;
	computed_likelihoods.clear();

	string line;
	size_t index = 0;
	while(getline(post_file, line)){
		vector<string> fields;
		string token;
		istringstream iss (line);
		while(getline(iss, token, '\t')) {
			fields.push_back(token);
		}

		REQUIRE(fields.size() == 6);

		// if variant has changed, compare likelihoods
		if ( (fields[0] != prev_var) && (index != 0)) {
			// compare computed likelihoods to expected ones
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 0));
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 1));
			computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 1, (unsigned char) 1));
			cleared = true;
			genotypes = GenotypingResult();
		}

		// determine which genotype the current state refers to
		unsigned char allele1 = (unsigned char) atoi(fields[4].c_str());
		unsigned char allele2 = (unsigned char) atoi(fields[5].c_str());
		genotypes.add_to_likelihood(allele1, allele2, string_to_l_double(fields[3]));
		cleared = false;
		index += 1;
		prev_var = fields[0];		
	
	}

	if (!cleared) {
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 0));
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 0, (unsigned char) 1));
		computed_likelihoods.push_back(genotypes.get_genotype_likelihood((unsigned char) 1, (unsigned char) 1));
	}

	// variants are enumerated in reverse order, that is why likelihoods also appear in reverse order (v2, v1)
	vector<double> expected_likelihoods_rev = { 0.9678020017, 0.031003181, 0.0011948172};

	REQUIRE(computed_likelihoods.size() == expected_likelihoods_rev.size());
	REQUIRE(compare_vectors(expected_likelihoods_rev, computed_likelihoods));
}