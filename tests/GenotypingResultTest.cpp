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

	// empty object, genotype should be unknown
	GenotypingResult r3;
	REQUIRE(r3.get_likeliest_genotype() == pair<int,int>(-1,-1));
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

TEST_CASE("GenotypingResult get_specific_likelihoods", "[GenotypingResult get_specific_likelihoods]") {
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.01);
	r.add_to_likelihood(0,1,0.02);
	r.add_to_likelihood(0,2,0.1);
	r.add_to_likelihood(1,0,0.15);
	r.add_to_likelihood(1,1,0.05);
	r.add_to_likelihood(1,2,0.15);
	r.add_to_likelihood(2,0,0.20);
	r.add_to_likelihood(2,1,0.22);
	r.add_to_likelihood(2,2,0.1);

	vector<long double> all_likelihoods = r.get_all_likelihoods(3);
	vector<long double> expected_likelihoods = {0.01, 0.17, 0.05, 0.30, 0.37, 0.1};

	for (size_t i = 0; i < 6; ++i) {
		REQUIRE(doubles_equal(all_likelihoods[i], expected_likelihoods[i]));
	}

	vector<unsigned char> defined_alleles = {0,2};
	GenotypingResult specific_likelihoods = r.get_specific_likelihoods(defined_alleles);
	vector<long double> likelihoods = specific_likelihoods.get_all_likelihoods(2);
	vector<long double> expected_specific_likelihoods = {0.0243902439, 0.73170731706, 0.24390243902};

	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(doubles_equal(likelihoods[i], expected_specific_likelihoods[i]));
	}
}

TEST_CASE("GenotypingResult get_specific_likelihoods2", "[GenotypingResult get_specific_likelihoods2]") {
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.2);
	r.add_to_likelihood(0,1,0.7);
	r.add_to_likelihood(1,1,0.1);

	vector<long double> all_likelihoods = r.get_all_likelihoods(2);
	vector<long double> expected_likelihoods = {0.2, 0.7, 0.1};

	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(doubles_equal(all_likelihoods[i], expected_likelihoods[i]));
	}

	vector<unsigned char> defined_alleles = {0,1};
	GenotypingResult specific_likelihoods = r.get_specific_likelihoods(defined_alleles);
	vector<long double> likelihoods = specific_likelihoods.get_all_likelihoods(2);

	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(doubles_equal(likelihoods[i], expected_likelihoods[i]));
	}
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

TEST_CASE("GenotypingResult combine", "[GenotypingResult combine]") {
	GenotypingResult r1;
	r1.add_to_likelihood(1,1,0.8);
	r1.add_to_likelihood(0,1,0.1);
	r1.add_to_likelihood(0,0,0.1);

	GenotypingResult r2;
	r2.add_to_likelihood(1,1,0.9);
	r2.add_to_likelihood(0,1,0.06);
	r2.add_to_likelihood(0,0,0.04);

	r1.combine(r2);
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(1,1), 1.7));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,1), 0.16));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,0), 0.14));

	GenotypingResult r3;
	r3.add_to_likelihood(1,2,0.1);
	r3.add_to_likelihood(0,1,0.2);
	r3.add_to_likelihood(0,0,0.4);
	r3.add_to_likelihood(2,2,0.3);

	r1.combine(r3);
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(2,2), 0.3));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(1,2), 0.1));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(1,1), 1.7));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,1), 0.36));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,0), 0.54));
	
}

TEST_CASE("GenotypingResult combine_empty1", "[GenotypingResult combine_empty1]") {
	// r1 is emty, all fields are 0
	GenotypingResult r1;
	GenotypingResult r2;
	r2.add_to_likelihood(1,1,0.9);
	r2.add_to_likelihood(0,1,0.06);
	r2.add_to_likelihood(0,0,0.04);

	r1.combine(r2);
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(1,1), 0.9));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,1), 0.06));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,0), 0.04));
}

TEST_CASE("GenotypingResult combine_empty2", "[GenotypingResult combine_empty2]") {
	// r2 is emty, all fields are 0
	GenotypingResult r1;
	r1.add_to_likelihood(1,1,0.9);
	r1.add_to_likelihood(0,1,0.06);
	r1.add_to_likelihood(0,0,0.04);
	GenotypingResult r2;

	r2.combine(r1);
	REQUIRE(doubles_equal(r2.get_genotype_likelihood(1,1), 0.9));
	REQUIRE(doubles_equal(r2.get_genotype_likelihood(0,1), 0.06));
	REQUIRE(doubles_equal(r2.get_genotype_likelihood(0,0), 0.04));
}

TEST_CASE("GenotypingResult combine_empty3", "[GenotypingResult combine_empty3]") {
	GenotypingResult r1;
	GenotypingResult r2;

	r1.combine(r2);
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(1,1), 0.0));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,1), 0.0));
	REQUIRE(doubles_equal(r1.get_genotype_likelihood(0,0), 0.0));
}

TEST_CASE("GenotypingResult normalize", "[GenotypingResult normalize]") {
	GenotypingResult g;
	g.add_to_likelihood(1,1,2);
	g.add_to_likelihood(1,0,1);
	g.add_to_likelihood(0,0,2);
	g.normalize();

	REQUIRE(doubles_equal(g.get_genotype_likelihood(1,1), 0.4));
	REQUIRE(doubles_equal(g.get_genotype_likelihood(0,1), 0.2));
	REQUIRE(doubles_equal(g.get_genotype_likelihood(0,0), 0.4));
}

TEST_CASE("GenotypingResult coverage_kmers", "[GenotypingResult coverage_kmers]") {
	GenotypingResult g;
	g.set_coverage(30);
	REQUIRE(g.coverage() == 30);
	g.set_unique_kmers(300);
	REQUIRE(g.nr_unique_kmers() == 300);
}
