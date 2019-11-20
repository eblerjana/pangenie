#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include "../src/collapsedhmm.hpp"
#include "../src/probabilitycomputer.hpp"
#include "utils.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("CollapsedHMM get_genotyping_result", "[CollapsedHMM get_genotyping_result]") {
	UniqueKmers u1(0,2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);

	UniqueKmers u2(1,3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(CopyNumber(0.01,0.01,0.9), a1);
	u2.insert_kmer(CopyNumber(0.9,0.3,0.1), a2);

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("NNN", "NNN", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,1}));
//	variants.push_back(Variant("NNN", "NNN", "chr1", 3000, 3003, {"CCC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	cout << "Before constructing HMM" << endl;
	CollapsedHMM hmm (&unique_kmers, true, true, 446.287102628, false, 0.25);
	cout << "After constructing HMM" << endl;

	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 2);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("CollapsedHMM get_genotyping_result_normalized", "[CollapsedHMM get_genotyping_result_normalized]") {
	UniqueKmers u1(0,2000);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1,0.0), a1);
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1,0.0), a2);

	UniqueKmers u2(1,3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(CopyNumber(0.01,0.01,0.9,0.0), a1);
	u2.insert_kmer(CopyNumber(0.9,0.3,0.1,0.0), a2);

	vector<UniqueKmers*> unique_kmers = {&u1,&u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("NNN", "NNN", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,1}));
//	variants.push_back(Variant("NNN", "NNN", "chr1", 3000, 3003, {"CCC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	CollapsedHMM hmm (&unique_kmers, true, true, 446.287102628, false, 0.25);
	
	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = { 0.0509465435, 0.9483202731, 0.0007331832, 0.9678020017, 0.031003181, 0.0011948172 };
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 2);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}


TEST_CASE("CollapsedHMM no_alt_allele", "[CollapsedHMM no_alt_allele]") {
	UniqueKmers u(0,2000);
	vector<unsigned char> a1 = {0,1};
	vector<unsigned char> a2 = {};
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,0);

	u.insert_kmer(CopyNumber(0.1,0.2,0.9), a1);
	u.insert_kmer(CopyNumber(0.3,0.4,0.1), a2);

	cout << "checkpoint 1" << endl;

	vector<UniqueKmers*> unique_kmers = {&u};
//	vector<Variant> variants;
//	variants.push_back(Variant("ATGC", "TGGG", "chr1", 2000, 2003, {"AAT", "ATT"}, {0,0,0}));
	CollapsedHMM hmm (&unique_kmers, true, true, 1.26, false, 0.25);

	cout << "checkpoint 2" << endl;

	// since only ref allele is covered by paths, 0/0 should have prob. 1
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,0), 1.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,1), 0.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(1,1), 0.0));
}

TEST_CASE("CollapsedHMM no_ref_allele", "[CollapsedHMM no_ref_allele]") {
	UniqueKmers u (0,2000);
	vector<unsigned char> a1 = {0,1};
	vector<unsigned char> a2 = {};
	u.insert_path(0,1);
	u.insert_path(1,1);
	u.insert_path(2,1);
	u.insert_kmer(CopyNumber(0.1,0.2,0.9), a1);
	u.insert_kmer(CopyNumber(0.3,0.4,0.1), a2);

	vector<UniqueKmers*> unique_kmers = {&u};
//	vector<Variant> variants;
//	variants.push_back(Variant("ATGC", "TGGG", "chr1", 2000, 2003, {"AAT", "ATT"}, {1,1,1}));
	CollapsedHMM hmm (&unique_kmers, true, true, 1.26, false, 0.25);
	// since only alt allele is covered by paths, 1/1 should have genotype 1/1
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(1,1), 1.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,1), 0.0));
	REQUIRE(doubles_equal(hmm.get_genotyping_result()[0].get_genotype_likelihood(0,0), 0.0));
}

TEST_CASE("CollapsedHMM no_unique_kmers", "[CollapsedHMM no_unique_kmers]") {
	UniqueKmers u1(0, 2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);

	UniqueKmers u2(0, 3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "TTT", "chr1", 2000, 2003, {"ATA", "AAA"}, {0,1}));
//	variants.push_back(Variant("CCC", "GGG", "chr1", 3000, 3003, {"CTC", "CGC"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	CollapsedHMM hmm (&unique_kmers, true, true, 446.287102628, false, 0.25);

	// each path combination should be equally likely here
	vector<double> expected_likelihoods = {0.25, 0.5, 0.25, 0.25, 0.5, 0.25};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 0);
	}

	REQUIRE( compare_vectors(computed_likelihoods, expected_likelihoods) );
}

TEST_CASE("CollapsedHMM no_unique_kmers2", "[CollapsedHMM no_unique_kmers2]") {
	UniqueKmers u1(0, 2000);
	u1.insert_path(0,0);
	u1.insert_path(1,0);
	u1.insert_path(2,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	UniqueKmers u2(0, 3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_path(2,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "TTT", "chr1", 2000, 2003, {"ATA", "AAA"}, {0,0,1}));
//	variants.push_back(Variant("CCC", "GGG", "chr1", 3000, 3003, {"CTC", "CGC"}, {0,1,1}));

	// recombination rate leads to recombination probability of 0.1
	CollapsedHMM hmm (&unique_kmers, true, true, 1070.02483182, false, 0.25);

	// each path combination should be equally likely here
	vector<double> expected_likelihoods = {4.0/9.0, 4.0/9.0, 1.0/9.0, 1.0/9.0, 4.0/9.0, 4.0/9.0};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 0);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("CollapsedHMM no_unique_kmers3", "[CollapsedHMM no_unique_kmers3]") {
	UniqueKmers u1(0,2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
	u1.insert_kmer(CopyNumber(0.1,0.9,0.1), a2);

	// no unique kmers for second variant
	UniqueKmers u2(1,3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	UniqueKmers u3(2,4000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
	u3.insert_kmer(CopyNumber(0.1,0.9,0.1), a1);
	u3.insert_kmer(CopyNumber(0.1,0.8,0.1), a2);

	vector<UniqueKmers*> unique_kmers = {&u1,&u2,&u3};
//	vector<Variant> variants;
//	variants.push_back(Variant("AAA", "CCC", "chr1", 2000, 2003, {"AGG", "ATG"}, {0,1}));
//	variants.push_back(Variant("GGG", "TCT", "chr1", 3000, 3003, {"TTT", "TGT"}, {0,1}));
//	variants.push_back(Variant("ATT", "TTC", "chr1", 4000, 4003, {"CAT", "CGT"}, {0,1}));

	// recombination rate leads to recombination probability of 0.1
	CollapsedHMM hmm (&unique_kmers, true, true, 446.287102628, false, 0.25);
	
	// expected likelihoods, as computed by hand
	vector<double> expected_likelihoods = {0.00264169937, 0.99471660125, 0.00264169937, 0.02552917716, 0.94894164567, 0.02552917716, 0.002961313333, 0.99407737333, 0.002961313333};
	vector<unsigned char> expected_haplotype1 = {0,0,0};
	vector<unsigned char> expected_haplotype2 = {1,1,1};
	vector<double> computed_likelihoods;
	vector<unsigned char> computed_haplotype1;
	vector<unsigned char> computed_haplotype2;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		pair<unsigned char,unsigned char> ht = result.get_haplotype();
		computed_haplotype1.push_back(ht.first);
		computed_haplotype2.push_back(ht.second);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods));

	// order of haplotype sequences can be different
	REQUIRE( ( ((expected_haplotype1 == computed_haplotype1) && (expected_haplotype2 == computed_haplotype2)) || ((expected_haplotype1 == computed_haplotype2) && (expected_haplotype2 == computed_haplotype1))) );
}

TEST_CASE("CollapsedHMM no_unique_kmers_uniform", "[CollapsedHMM no_unique_kmers_uniorm]") {
	UniqueKmers u1 (0,2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_path(2,1);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);

	UniqueKmers u2 (0,3000);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_path(2,1);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2};
	// switch on uniform transition probabilities
	CollapsedHMM hmm (&unique_kmers, true, true, 1.26, true, 0.25);

	vector<double> expected_likelihoods = {1/9.0, 4/9.0, 4/9.0, 4/9.0, 4/9.0, 1/9.0};
	vector<double> computed_likelihoods;
	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 0);
	}

	for (size_t i = 0; i < expected_likelihoods.size(); ++i) {
		cout << expected_likelihoods[i] << " " << computed_likelihoods[i] << endl;
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("CollapsedHMM only_kmers", "[CollapsedHMM only_kmers]") {
	/** PGGTyper-kmers: insert one ref and one alt path and allow uniform transitions between paths **/
	UniqueKmers u1 (0,2000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(CopyNumber(0.05,0.9,0.05), a1);
	u1.insert_kmer(CopyNumber(0.1,0.7,0.2), a2);

	UniqueKmers u2 (1,3000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(CopyNumber(0.9,0.07,0.03), a1);
	u2.insert_kmer(CopyNumber(0.1,0.2,0.7), a2);

	UniqueKmers u3 (2, 4000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
	u3.insert_kmer(CopyNumber(0.6,0.3,0.1), a1);
	u3.insert_kmer(CopyNumber(0.3,0.4,0.3), a2);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	// switch on uniform transition probabilities
	CollapsedHMM hmm (&unique_kmers, true, true, 1.26, true, 0.25);

	vector<double> expected_likelihoods = {0.00392156862745098, 0.988235294117647, 0.00784313725490196, 0.0045385779122541605, 0.02118003025718608, 0.9531013615733737, 0.06666666666666667, 0.5333333333333333, 0.39999999999999997};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 2);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("CollapsedHMM emissions_zero", "[CollapsedHMM emissions_zero]") {
	UniqueKmers u1 (0,1000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);

	UniqueKmers u2(0,2000);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a1);
	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a1);

	UniqueKmers u3(0,3000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	CollapsedHMM hmm (&unique_kmers, true, true, 446.287102628, false, 0.25);
	vector<double> expected_likelihoods = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 2);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}

TEST_CASE("CollapsedHMM underflow", "[CollapsedHMM underflow]") {
	UniqueKmers u1 (0,1000);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
	u1.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);

	UniqueKmers u2 (0,2000);
	u2.insert_path(0,0);
	u2.insert_path(1,1);
	u2.insert_kmer(CopyNumber(0.0,0.0,1.0), a1);
	u2.insert_kmer(CopyNumber(1.0,0.0,0.0), a2);

	UniqueKmers u3 (0,3000);
	u3.insert_path(0,0);
	u3.insert_path(1,1);
	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a1);
	u3.insert_kmer(CopyNumber(0.0,1.0,0.0), a2);

	vector<UniqueKmers*> unique_kmers = {&u1, &u2, &u3};
	CollapsedHMM hmm (&unique_kmers, true, true, 0.0, false, 0.25);
	vector<double> expected_likelihoods = {0.0,1.0,0.0,0.25,0.5,0.25,0.0,1.0,0.0};
	vector<double> computed_likelihoods;

	for (auto result : hmm.get_genotyping_result()) {
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,0));
		computed_likelihoods.push_back(result.get_genotype_likelihood(0,1));
		computed_likelihoods.push_back(result.get_genotype_likelihood(1,1));
		REQUIRE(result.get_nr_unique_kmers() == 2);
	}

	REQUIRE( compare_vectors(expected_likelihoods, computed_likelihoods) );
}
