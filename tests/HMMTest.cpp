#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include "../src/hmm.hpp"
#include "../src/probabilitycomputer.hpp"
#include "utils.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("HMM get_genotyping_result", "[HMM get_genotyping_result]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);
	u1.set_coverage(5);

	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(5, a2);
	u2.set_coverage(5);

	ProbabilityTable probs(5, 10, 30, 0.0L);
	probs.modify_probability(5, 10, CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(5, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(5, 5, CopyNumber(0.9,0.3,0.1));
	vector<UniqueKmers*> unique_kmers = {&u1,&u2};

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == 1);
//		REQUIRE(result.get_allele_kmer_count(1) == 1);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("HMM skip_reference_position", "[HMM skip_reference_position]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);
	u1.set_coverage(5);

	// position only covered by reference alleles and should be skipped
	UniqueKmers u2(2500);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_kmer(10, a1);
	u2.insert_kmer(20, a2);

	UniqueKmers u3(3000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
	u3.insert_kmer(20, a1);
	u3.insert_kmer(5, a2);
	u3.set_coverage(5);

	ProbabilityTable probs(5, 10, 30, 0.0L);
	probs.modify_probability(5, 10, CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(5, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(5, 5, CopyNumber(0.9,0.3,0.1));
	vector<UniqueKmers*> unique_kmers = {&u1,&u2,&u3};

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.0, 0.0, 0.0, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	unsigned int index = 0;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		if (index != 1) {
//			REQUIRE(result.get_nr_unique_kmers() == 2);
//			REQUIRE(result.get_allele_kmer_count(0) == 1);
//			REQUIRE(result.get_allele_kmer_count(1) == 1);
		} else {
//			REQUIRE(result.get_nr_unique_kmers() == 0);
//			REQUIRE(result.get_allele_kmer_count(0) == 0);
//			REQUIRE(result.get_allele_kmer_count(1) == 0);
		}
		index += 1;
	}
	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}


TEST_CASE("HMM get_genotyping_result_normalized", "[HMM get_genotyping_result_normalized]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(1, a2);
	
	ProbabilityTable probs (0, 5, 30, 0.0L);
	probs.modify_probability(0, 10, CopyNumber(0.1,0.9,0.1,0.0));
	probs.modify_probability(0, 20, CopyNumber(0.01,0.01,0.9,0.0));
	probs.modify_probability(0, 1, CopyNumber(0.9,0.3,0.1,0.0));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);
	
	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == 1);
//		REQUIRE(result.get_allele_kmer_count(1) == 1);
	}
	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("HMM undefined_alleles1", "[HMM get_undefined_alleles1]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_empty_allele(0, true);
	u1.insert_empty_allele(1);
	REQUIRE (u1.is_undefined_allele(0));
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(10, a1);

	UniqueKmers u2(3000);
	u2.insert_path(0,1);
	u2.insert_path(1,0);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(1, a2);

	ProbabilityTable probs(0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(0,20,CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0,1,CopyNumber(0.9,0.3,0.1));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
	// recombination rate leads to recombination probability of 0.1
	// TODO: this currently runs ONLY the genotyping
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.02396597038, 0.52185641164,  0.45417761795, 0.97855858361, 0.01875778106,  0.00268363531};
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );

	computed_likelihoods.clear();
	// extract only likelihoods for defined genotypes
	vector<double> expected_specific_likelihoods = {1.0, 0.0, 0.0, 0.97855858361, 0.01875778106, 0.00268363531};
	vector<vector<unsigned char>> defined_alleles = {{1}, {0,1}};
	vector<GenotypingResult> result = hmm.get_genotyping_result();
	for (unsigned int i = 0; i < result.size(); ++i) {
		GenotypingResult final_likelihood = result[i].get_specific_likelihoods(defined_alleles[i]);
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(1,1));
	}

	REQUIRE( compare_vectors(expected_specific_likelihoods, computed_likelihoods) );
}

TEST_CASE("HMM undefined_alleles2", "[HMM get_undefined_alleles1]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_path(0,0);
	u1.insert_path(1,0);

	UniqueKmers u2(3000);
	u2.insert_empty_allele(0, true);
	u2.insert_empty_allele(1);
	REQUIRE (u2.is_undefined_allele(0));
	u2.insert_path(0,1);
	u2.insert_path(1,0);
	u2.insert_kmer(20, a2);
	u2.insert_kmer(1, a1);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0, 1, CopyNumber(0.9,0.3,0.1));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
	// recombination rate leads to recombination probability of 0.1
	// TODO: this currently runs ONLY the genotyping
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = {0.0,0.0,0.0,0.11813512445,0.1617937574,0.72007111814};
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );

	computed_likelihoods.clear();
	// extract only likelihoods for defined genotypes
	vector<double> expected_specific_likelihoods = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
	vector<vector<unsigned char>> defined_alleles = {{0,1}, {0}};
	vector<GenotypingResult> result = hmm.get_genotyping_result();
	for (unsigned int i = 0; i < result.size(); ++i) {
		GenotypingResult final_likelihood = result[i].get_specific_likelihoods(defined_alleles[i]);
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(1,1));
	}

	REQUIRE( compare_vectors(expected_specific_likelihoods, computed_likelihoods) );

}

TEST_CASE("HMM only_undefined_alleles", "[HMM only_undefined_alleles]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	// set both alleles to undefined
	u1.insert_empty_allele(0, true);
	u1.insert_empty_allele(1, true);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	UniqueKmers u2(3000);
	// set both alleles to undefined
	u2.insert_empty_allele(0, true);
	u2.insert_empty_allele(1, true);
	u2.insert_path(0,1);
	u2.insert_path(1,0);
//	u2.insert_kmer(CopyNumber(0.01,0.01,0.9), a1);
//	u2.insert_kmer(CopyNumber(0.9,0.3,0.1), a2);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(1, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(0,20,CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0,1,CopyNumber(0.9,0.3,0.1));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
	// recombination rate leads to recombination probability of 0.1
	// TODO: this currently runs ONLY the genotyping
	HMM hmm (&unique_kmers, &probs, true, false, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 0);
//		REQUIRE(result.get_allele_kmer_count(0) == 0);
//		REQUIRE(result.get_allele_kmer_count(1) == 0);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );

	computed_likelihoods.clear();
	// extract only likelihoods for defined genotypes
	vector<double> expected_specific_likelihoods = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	vector<vector<unsigned char>> defined_alleles = {{}, {}};
	vector<GenotypingResult> result = hmm.get_genotyping_result();
	for (unsigned int i = 0; i < result.size(); ++i) {
		GenotypingResult final_likelihood = result[i].get_specific_likelihoods(defined_alleles[i]);
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(final_likelihood.get_genotype_likelihood(1,1));
	}

	REQUIRE( compare_vectors(expected_specific_likelihoods, computed_likelihoods) );
}


TEST_CASE("HMM no_alt_allele", "[HMM no_alt_allele]") {
	UniqueKmers u(2000);
	vector<unsigned char> a1 = {0,1};
	vector<unsigned char> a2 = {};
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,0);

//	u.insert_kmer(CopyNumber(0.1,0.2,0.9), a1);
//	u.insert_kmer(CopyNumber(0.3,0.4,0.1), a2);
	u.insert_kmer(10, a1);
	u.insert_kmer(5, a2);

	ProbabilityTable probs (0,1,11, 0.0L);
	probs.modify_probability(0,10,CopyNumber(0.1,0.2,0.9));
	probs.modify_probability(0,5,CopyNumber(0.3,0.4,0.1));

	vector<UniqueKmers*> unique_kmers = {&u};
//	vector<Variant> variants;
//	variants.push_back(Variant("ATGC", "TGGG", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,0,0}));
	HMM hmm (&unique_kmers, &probs, true, true, 1.26, false, 0.25);
	// since only ref allele is covered by paths, all likelihoods are set to zero (column is not considered).
	REQUIRE(hmm.get_genotyping_result()[0].get_likeliest_genotype() == pair<int,int>(-1,-1));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,0), 0.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,1), 0.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(1,1), 0.0));
//	REQUIRE(hmm.get_genotyping_result()[0].get_allele_kmer_count(0) == 0);
//	REQUIRE(hmm.get_genotyping_result()[0].get_allele_kmer_count(1) == 0);
}



TEST_CASE("HMM no_ref_allele", "[HMM no_ref_allele]") {
	UniqueKmers u (2000);
	vector<unsigned char> a1 = {0,1};
	vector<unsigned char> a2 = {};
	u.insert_path(0,1);
	u.insert_path(1,1);
	u.insert_path(2,1);
//	u.insert_kmer(CopyNumber(0.1,0.2,0.9), a1);
//	u.insert_kmer(CopyNumber(0.3,0.4,0.1), a2);
	u.insert_kmer (20, a1);
	u.insert_kmer (10, a2);

	ProbabilityTable probs(0, 1, 21, 0.0L);
	probs.modify_probability(0, 20, CopyNumber(0.1,0.2,0.9));
	probs.modify_probability(0, 10, CopyNumber(0.3,0.4,0.1));

	vector<UniqueKmers*> unique_kmers = {&u};
//	vector<Variant> variants;
//	variants.push_back(Variant("ATGC", "TGGG", "chr1", 2000, 2003, {"AAT", "ATT"}, {1,1,1}));
	HMM hmm (&unique_kmers, &probs, true, true, 1.26, false, 0.25);
	// since only alt allele is covered by paths, 1/1 should have genotype 1/1
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(1,1), 1.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,1), 0.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,0), 0.0));
//	REQUIRE(hmm.get_genotyping_result()[0].get_allele_kmer_count(0) == 1);
//	REQUIRE(hmm.get_genotyping_result()[0].get_allele_kmer_count(1) == 1);
}


TEST_CASE("HMM no_unique_kmers", "[HMM no_unique_kmers]") {
	UniqueKmers u1(2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);

	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	ProbabilityTable probs;

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "TTT", "chr1", 2000, 2003, {"ATA", "AAA"}, {0,1}));
//	variants.push_back(Variant("CCC", "GGG", "chr1", 3000, 3003, {"CTC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// each path combination should be equally likely here
	vector<double> expected_likelihoods = {0.25, 0.5, 0.25, 0.25, 0.5, 0.25};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 0);
//		REQUIRE(result.get_allele_kmer_count(0) == 0);
//		REQUIRE(result.get_allele_kmer_count(1) == 0);
	}

	REQUIRE( compare_vectors(computed_likelihoods, expected_likelihoods) );
}



TEST_CASE("HMM no_unique_kmers2", "[HMM no_unique_kmers2]") {
	UniqueKmers u1(2000);
	u1.insert_path(0,0);
	u1.insert_path(1,0);
	u1.insert_path(2,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_path(2,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	ProbabilityTable probs;

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "TTT", "chr1", 2000, 2003, {"ATA", "AAA"}, {0,0,1}));
//	variants.push_back(Variant("CCC", "GGG", "chr1", 3000, 3003, {"CTC", "CGC"}, {0,1,1}));

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 1070.02483182, false, 0.25);

	// each path combination should be equally likely here
	vector<double> expected_likelihoods = {4.0/9.0, 4.0/9.0, 1.0/9.0, 1.0/9.0, 4.0/9.0, 4.0/9.0};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 0);
//		REQUIRE(result.get_allele_kmer_count(0) == 0);
//		REQUIRE(result.get_allele_kmer_count(1) == 0);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}



TEST_CASE("HMM no_unique_kmers3", "[HMM no_unique_kmers3]") {
	UniqueKmers u1(2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	// no unique kmers for second variant
	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	UniqueKmers u3(4000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
//	u3.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
//	u3.insert_kmer(CopyNumber(0.1,0.8,0.1), a2);
	u3.insert_kmer(10, a1);
	u3.insert_kmer(9, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0, 10, CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(0, 9, CopyNumber(0.1,0.8,0.1));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2,&u3};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "CCC", "chr1", 2000, 2003, {"AGG", "ATG"}, {0,1}));
//	variants.push_back(Variant("GGG", "TCT", "chr1", 3000, 3003, {"TTT", "TGT"}, {0,1}));
//	variants.push_back(Variant("ATT", "TTC", "chr1", 4000, 4003, {"CAT", "CGT"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);
	
	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = {0.00264169937, 0.99471660125, 0.00264169937, 0.02552917716, 0.94894164567, 0.02552917716, 0.002961313333, 0.99407737333, 0.002961313333};
	vector<unsigned char> expected_haplotype1 = {0,0,0};
	vector<unsigned char> expected_haplotype2 = {1,1,1};
	vector<vector<unsigned int>> expected_counts = { {1,1}, {0,0}, {1,1} };
	vector<double> computed_likelihoods;
	vector<unsigned char> computed_haplotype1;
	vector<unsigned char> computed_haplotype2;
	unsigned int index = 0;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		pair<unsigned char,unsigned char> ht = result.get_haplotype();
		computed_haplotype1.push_back(ht.first);
		computed_haplotype2.push_back(ht.second);
//		REQUIRE(result.get_allele_kmer_count(0) == expected_counts[index][0]);
//		REQUIRE(result.get_allele_kmer_count(1) == expected_counts[index][1]);
		index += 1;
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods));

	// order of haplotype sequences can be different
	REQUIRE( ( ((expected_haplotype1 == computed_haplotype1) && (expected_haplotype2 == computed_haplotype2)) || ((expected_haplotype1 == computed_haplotype2) && (expected_haplotype2 == computed_haplotype1))) );
}


TEST_CASE("HMM no_unique_kmers_uniform", "[HMM no_unique_kmers_uniorm]") {
	UniqueKmers u1 (2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_path(2,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);

	UniqueKmers u2 (3000);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_path(2,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	ProbabilityTable probs;

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
	// switch on uniform transition probabilities
	HMM hmm (&unique_kmers, &probs, true, true, 1.26, true, 0.25);

	vector<double> expected_likelihoods = {1/9.0, 4/9.0, 4/9.0, 4/9.0, 4/9.0, 1/9.0};
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 0);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}



TEST_CASE("HMM only_kmers", "[HMM only_kmers]") {
	// PGGTyper-kmers: insert one ref and one alt path and allow uniform transitions between paths
	UniqueKmers u1 (2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
//	u1.insert_kmer(CopyNumber(0.05,0.9,0.05), a1);
//	u1.insert_kmer(CopyNumber(0.1,0.7,0.2), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(12, a2);

	UniqueKmers u2 (3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
//	u2.insert_kmer(CopyNumber(0.9,0.07,0.03), a1);
//	u2.insert_kmer(CopyNumber(0.1,0.2,0.7), a2);
	u2.insert_kmer(1, a1);
	u2.insert_kmer(20, a2);

	UniqueKmers u3 (4000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
//	u3.insert_kmer(CopyNumber(0.6,0.3,0.1), a1);
//	u3.insert_kmer(CopyNumber(0.3,0.4,0.3), a2);
	u3.insert_kmer(5, a1);
	u3.insert_kmer(7, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.05,0.9,0.05));
	probs.modify_probability(0,12,CopyNumber(0.1,0.7,0.2));
	probs.modify_probability(0,1,CopyNumber(0.9,0.07,0.03));
	probs.modify_probability(0,20,CopyNumber(0.1,0.2,0.7));
	probs.modify_probability(0,5,CopyNumber(0.6,0.3,0.1));
	probs.modify_probability(0,7,CopyNumber(0.3,0.4,0.3));

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	// switch on uniform transition probabilities
	HMM hmm (&unique_kmers, &probs, true, true, 1.26, true, 0.25);

	vector<double> expected_likelihoods = {0.00392156862745098, 0.988235294117647, 0.00784313725490196, 0.0045385779122541605, 0.0423600605143722, 0.9531013615733737, 0.06666666666666667, 0.5333333333333333, 0.39999999999999997};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == 1);
//		REQUIRE(result.get_allele_kmer_count(1) == 1);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}



TEST_CASE("HMM emissions_zero", "[HMM emissions_zero]") {
	UniqueKmers u1 (1000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
//	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
//	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	UniqueKmers u2(2000);
	u2.insert_path(0,1);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
//	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a1);
//	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a1);
	u2.insert_kmer(0, a2);
	u2.insert_kmer(0, a2);

	UniqueKmers u3(3000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
//	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
//	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);
	u3.insert_kmer(10, a1);
	u3.insert_kmer(10, a2);

	ProbabilityTable probs (0,1,11,0.0L);
	probs.modify_probability(0, 10, CopyNumber(0.0,1.0,0.0));
	probs.modify_probability(0, 0, CopyNumber(1.0,0.0,0.0));

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);
	// currently, backward probabilities that all become zero, are not corrected in current column (but stored as uniform for further computations)
	// since all backward probs in first column are zero, likelihoods are zero as well. Backward probs are later set to uniform.
	vector<double> expected_likelihoods = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0};
	vector<double> computed_likelihoods;

	unsigned int index = 0;
	// second column is skipped
	vector<vector<unsigned int>> expected_counts = { {1,1}, {0,2}, {1,1} };

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == expected_counts[index][0]);
//		REQUIRE(result.get_allele_kmer_count(1) == expected_counts[index][1]); 
		index += 1;
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}


TEST_CASE("HMM underflow", "[HMM underflow]") {
	UniqueKmers u1 (1000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
//	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
//	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	UniqueKmers u2 (2000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
//	u2.insert_kmer(CopyNumber(0.0,0.0,1.0), a1);
//	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a2);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(0, a2);

	UniqueKmers u3 (3000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
//	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
//	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);
	u3.insert_kmer(10, a1);
	u3.insert_kmer(10, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.0,1.0,0.0));
	probs.modify_probability(0,20,CopyNumber(0.0,0.0,1.0));
	probs.modify_probability(0,0,CopyNumber(1.0,0.0,0.0));

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	HMM hmm (&unique_kmers, &probs, true, true, 0.0, false, 0.25);
	// currently, backward probabilities that all become zero, are not corrected in current column (but stored as uniform for further computations)
	vector<double> expected_likelihoods = {0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == 1);
//		REQUIRE(result.get_allele_kmer_count(1) == 1);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}


TEST_CASE("HMM get_genotyping_result_neutral_kmers", "[HMM get_genotyping_result_with_kmer]") {
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	vector<unsigned char> a3 = {0,1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);
	// these kmers are located on all alleles
//	u1.insert_kmer(CopyNumber(0.05, 0.45, 0.5), a3);
//	u1.insert_kmer(CopyNumber(0.4, 0.5, 0.1), a3);

	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);
	u1.insert_kmer(12, a3);
	u1.insert_kmer(5, a3);

	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
//	u2.insert_kmer(CopyNumber(0.01,0.01,0.9), a1);
//	u2.insert_kmer(CopyNumber(0.9,0.3,0.1), a2);
//	u2.insert_kmer(CopyNumber(0.01, 0.49, 0.5), a3);
//	u2.insert_kmer(CopyNumber(0.3, 0.4, 0.3), a3);

	u2.insert_kmer(20, a1);
	u2.insert_kmer(1, a2);
	u2.insert_kmer(15, a3);
	u2.insert_kmer(9, a3);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(0,12,CopyNumber(0.05, 0.45, 0.5));
	probs.modify_probability(0,5,CopyNumber(0.4, 0.5, 0.1));
	probs.modify_probability(0,20,CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0,1,CopyNumber(0.9,0.3,0.1));
	probs.modify_probability(0,15,CopyNumber(0.01, 0.49, 0.5));
	probs.modify_probability(0,9,CopyNumber(0.3, 0.4, 0.3));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("NNN", "NNN", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,1}));
//	variants.push_back(Variant("NNN", "NNN", "chr1", 3000, 3003, {"CCC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 4);
//		REQUIRE(result.get_allele_kmer_count(0) == 3);
//		REQUIRE(result.get_allele_kmer_count(1) == 3);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}



TEST_CASE("HMM only_paths", "[HMM only_paths]") {
	// want to only consider the following paths. All others should be ignored.
	vector<unsigned short> only_paths = {0,3};
	UniqueKmers u1(2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	vector<unsigned char> a3 = {2};
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_empty_allele(2);
	u1.insert_path(0,0);
	u1.insert_path(1,2);
	u1.insert_path(2,1);
	u1.insert_path(3,1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
//	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);
	u1.insert_kmer(10, a1);
	u1.insert_kmer(10, a2);

	UniqueKmers u2(3000);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);
	u2.insert_empty_allele(2);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_path(2,2);
	u2.insert_path(3,1);
//	u2.insert_kmer(CopyNumber(0.01,0.01,0.9), a1);
//	u2.insert_kmer(CopyNumber(0.9,0.3,0.1), a2);
	u2.insert_kmer(20, a1);
	u2.insert_kmer(1, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0,10,CopyNumber(0.1,0.9,0.1));
	probs.modify_probability(0,20,CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0,1,CopyNumber(0.9,0.3,0.1));

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("NNN", "NNN", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,1}));
//	variants.push_back(Variant("NNN", "NNN", "chr1", 3000, 3003, {"CCC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25, &only_paths);

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 2);
//		REQUIRE(result.get_allele_kmer_count(0) == 1);
//		REQUIRE(result.get_allele_kmer_count(1) == 1);
//		REQUIRE(result.get_allele_kmer_count(2) == 0);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("HMM no_only_paths2", "[HMM only_paths2]") {
	vector<unsigned char> a2 = {2};
	UniqueKmers u1(2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_path(2,2);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_empty_allele(2);
//	u1.insert_kmer(CopyNumber(0.05, 0.8, 0.15), a2);
	u1.insert_kmer(12, a2);

	UniqueKmers u2(3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);
	u2.insert_empty_allele(2);

//	u2.insert_kmer(CopyNumber(0.05, 0.8, 0.15), a2);
	u2.insert_kmer(12,a2);
	u2.insert_path(2,2);

	ProbabilityTable probs (0,1,13,0.0L);
	probs.modify_probability(0,12,CopyNumber(0.05, 0.8, 0.15));

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "TTT", "chr1", 2000, 2003, {"ATA", "AAA"}, {0,1}));
//	variants.push_back(Variant("CCC", "GGG", "chr1", 3000, 3003, {"CTC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	vector<unsigned short> only_paths = {0,1};
	HMM hmm (&unique_kmers, &probs, true, true, 446.287102628, false, 0.25, &only_paths);

	// each path combination should be equally likely here
	vector<double> expected_likelihoods = {0.25, 0.5, 0.25, 0.25, 0.5, 0.25};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
//		REQUIRE(result.get_nr_unique_kmers() == 1);
//		REQUIRE(result.get_allele_kmer_count(0) == 0);
//		REQUIRE(result.get_allele_kmer_count(1) == 0);
	}

	REQUIRE( compare_vectors(computed_likelihoods, expected_likelihoods) );
}
