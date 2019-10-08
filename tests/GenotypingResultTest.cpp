#include "catch.hpp"
#include "utils.hpp"
#include "../src/genotypingresult.hpp"
#include <vector>
#include <string>

#include <iostream>

using namespace std;

TEST_CASE("GenotypingResult get_genotype_likelihood", "[GenotypingResult get_genotype_likelihood]"){
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.1);
	r.add_to_likelihood(0,0,0.2);
	r.add_to_likelihood(0,1,0.1);
	r.add_to_likelihood(1,0,0.1);
	r.add_to_likelihood(1,1,0.5);

	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,0), 0.3));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,1), 0.2));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(1,1), 0.5));

	r.add_first_haplotype_allele(0);
	r.add_second_haplotype_allele(1);

	REQUIRE(r.get_haplotype() == pair<unsigned char, unsigned char>(0,1));
	REQUIRE(r.get_likeliest_genotype() == pair<int, int>(1,1));
}

TEST_CASE("GenotypingResult get_likeliest_genotype", "[GenotypingResult get_likeliest_genotype]") {
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.1);
	r.add_to_likelihood(0,1,0.5);
	r.add_to_likelihood(1,1,0.4);
	REQUIRE(r.get_likeliest_genotype() == pair<int, int>(0,1));

	GenotypingResult r2;
	r2.add_to_likelihood(0,0,0.5);
	r2.add_to_likelihood(0,1,0.5);
	// both genotypes equally likely
	REQUIRE(r2.get_likeliest_genotype() == pair<int, int>(-1,-1));
}

TEST_CASE("GenotypingResult divide_likelihoods_by", "[GenotypingResult divide_likelihoods_by]") {
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.2);
	r.add_to_likelihood(0,1,0.8);
	r.add_to_likelihood(1,1,1.0);

	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,0), 0.2));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,1), 0.8));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(1,0), 0.8));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(1,1), 1.0));

	r.divide_likelihoods_by(2.0);
	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,0), 0.1));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(0,1), 0.4));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(1,0), 0.4));
	REQUIRE(doubles_equal(r.get_genotype_likelihood(1,1), 0.5));
}

TEST_CASE("GenotypingResult get_all_likelihoods (biallelic)", "[GenotypingResult get_all_likelihoods (biallelic)]") {
	
	// biallelic site
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.1);
	r.add_to_likelihood(1,1,0.2);
	r.add_to_likelihood(0,1,0.7);

	vector<long double> computed_likelihoods = r.get_all_likelihoods(2);
	vector<long double> expected_likelihoods = {0.1,0.7,0.2};
	REQUIRE(computed_likelihoods.size() == 3);

	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(doubles_equal(computed_likelihoods[i], expected_likelihoods[i]));
	}
}

TEST_CASE("GenotypingResult get_all_likelihoods (triallelic)", "[GenotypingResult get_all_likelihoods (triallelic]") {

	// triallelic site
	GenotypingResult r;
	r.add_to_likelihood(0,1,0.01);
	r.add_to_likelihood(0,0,0.05);
	r.add_to_likelihood(1,1,0.04);
	r.add_to_likelihood(2,2,0.3);
	r.add_to_likelihood(1,2,0.5);
	r.add_to_likelihood(0,2,0.1);

	vector<long double> computed_likelihoods = r.get_all_likelihoods(3);
	vector<long double> expected_likelihoods = {0.05,0.01,0.04,0.1,0.5,0.3};
	REQUIRE(computed_likelihoods.size() == 6);

	for (size_t i = 0; i < 6; ++i) {
		REQUIRE(doubles_equal(computed_likelihoods[i], expected_likelihoods[i]));
	}

	REQUIRE(r.get_genotype_quality(1,2) == 3);
	REQUIRE(r.get_genotype_quality(0,1) == 0);
}

TEST_CASE("GenotypingResult get_genotype_quality unnormalized", "[GenotypingResult get_genotype_quality unnormalized]") {

	GenotypingResult r;
	r.add_to_likelihood(0,0,0.4);
	r.add_to_likelihood(0,1,0.6);
	r.add_to_likelihood(1,1,0.7);

	CHECK_THROWS(r.get_genotype_quality(1,1));
	
	r.divide_likelihoods_by(1.7);
	REQUIRE(r.get_genotype_quality(1,1) == 2);
}

TEST_CASE("GenotypingResult get_genotype_quality", "[GenotypingResult get_genotype_quality]") {
	GenotypingResult r;
	r.add_to_likelihood(1,1,1.0);
	REQUIRE(r.get_genotype_quality(1,1) == 10000);
}

TEST_CASE("GenotypingResult get_nr_unique_kmers", "[GenotypingResult get_nr_unique_kmers]") {
	GenotypingResult r;
	r.set_nr_unique_kmers(4);
	REQUIRE(r.get_nr_unique_kmers() == 4);
	r.set_nr_unique_kmers(5);
	REQUIRE(r.get_nr_unique_kmers() == 5);
}
