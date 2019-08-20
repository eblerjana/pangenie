#include "catch.hpp"
#include "../src/variant.hpp"
#include <vector>
#include <string>
#include "utils.hpp"

using namespace std;

TEST_CASE("Variant testcase 1", "[Variant testcase 1]"){
	Variant v1("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATT"}, {0,1});
	Variant v2("GCT", "CCC", "chr1", 15, 16, {"A", "G"}, {1,0});

	REQUIRE(v1.nr_of_alleles() == 2);
	REQUIRE(v2.nr_of_alleles() == 2);
	REQUIRE(v1.nr_of_paths() == 2);
	REQUIRE(v2.nr_of_paths() == 2);
	REQUIRE(v1.get_allele_string(0) == "ATGC");
	REQUIRE(v2.get_allele_string(1) == "G");

	REQUIRE(v1.allele_on_path(0,0));
	REQUIRE(! v1.allele_on_path(1,0));
	REQUIRE(v2.allele_on_path(1,0));
	REQUIRE(v2.allele_on_path(0,1));

	v1.combine_variants(v2);
	REQUIRE(v1.get_allele_string(0) == "ATGCTA");
	REQUIRE(v1.get_allele_string(1) == "ATGCTG");
	REQUIRE(v1.get_allele_string(2) == "ATTTA");
	REQUIRE(v1.nr_of_alleles() == 3);
	REQUIRE(v1.is_combined());
	
	v1.add_flanking_sequence();
	REQUIRE(v1.get_allele_string(0) == "AAAATGCTACCC");
	REQUIRE(v1.get_allele_string(1) == "AAAATGCTGCCC");
	REQUIRE(v1.get_allele_string(2) == "AAAATTTACCC");

	v1.remove_flanking_sequence();
	REQUIRE(v1.get_allele_string(0) == "ATGCTA");
	REQUIRE(v1.get_allele_string(1) == "ATGCTG");
	REQUIRE(v1.get_allele_string(2) == "ATTTA");
}

TEST_CASE("Variant operator==", "[Variant operator==]") {
	Variant v1("AAA", "TAC", "chr1", 10, 13, {"ATG", "C"}, {0,1});
	Variant v2("GCT", "TTT", "chr1", 10, 13, {"ATG", "C"}, {0,1});

	// different flanks
	REQUIRE(!(v1 == v2));
	REQUIRE(v1 != v2);

	Variant v3("AAA", "TAC", "chr1", 10, 13, {"ATG", "CG"}, {0,1});
	// different alleles
	REQUIRE(!(v1 == v3));
	REQUIRE(v1 != v3);

	Variant v4("AAA", "TAC", "chr2", 10, 13, {"ATG", "CG"}, {0,1});
	// different chromosome
	REQUIRE(!(v1 == v4));
	REQUIRE(v1 != v4);

	Variant v5("AAA", "TAC", "chr2", 10, 12, {"AT", "CG"}, {0,1});
	// different positions
	REQUIRE(!(v4 == v5));
	REQUIRE(v4 != v5);

	Variant v6("AAA", "TAC", "chr1", 10, 13, {"ATG", "C"}, {1,0});
	// different paths
	REQUIRE(!(v2 == v6));
	REQUIRE(v2 != v6);

	Variant v7("AAA", "TAC", "chr1", 10, 13, {"ATG", "C"}, {0,1});
	REQUIRE( v1 == v7);
	REQUIRE(! (v1 != v7));

	v7.add_flanking_sequence();
	// flanks added
	REQUIRE(! (v1 == v7));
	REQUIRE(v1 != v7);
}


TEST_CASE("Variant invalid1", "[Variant invalid1]"){
	// flanking sequences of different sizes are not allowed
	CHECK_THROWS(Variant ("AAA", "TTAA", "chr1", 10, 14, {"ATGC", "ATT"}, {0,1}));
	// end position does not match reference allele's length
	CHECK_THROWS(Variant ("AAA", "TAA", "chr1", 10, 11, {"ATGC", "ATT"}, {0,1}));
	// path allele index does not match number of alleles
	CHECK_THROWS(Variant ("AAA", "TAA", "chr1", 10, 11, {"ATGC", "ATT"}, {0,2}));
	// start position larger than end position
	CHECK_THROWS(Variant ("AAA", "TAA", "chr1", 14, 10, {"ATGC", "ATT"}, {0,1}));
}


TEST_CASE("Variant combine_variants_invalid", "[Variant combine_variants_invalid]"){
	Variant v1("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATT"}, {0,1});
	Variant v2("TGCT", "CCCC", "chr1", 15, 16, {"A", "G"}, {1,0});
	// flanking sequences of different size
	CHECK_THROWS(v1.combine_variants(v2));
	REQUIRE(!v1.is_combined());
	Variant v3("CCC", "TTT", "chr2", 17, 18, {"A","G"}, {1,0});
	// variants too far apart
	CHECK_THROWS(v2.combine_variants(v3));
	REQUIRE(!v2.is_combined());
}

TEST_CASE("Variant combine_variants", "[Variant combine_variants]") {
	Variant v1 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	Variant v2 ("AACT", "ACTG", "chr2", 7, 10, {"GAG", "ACC"}, {0,0,1,1});
	Variant v3 ("GACT", "GGAA", "chr2", 13, 14, {"G", "GTC"}, {0,0,1,0});

	v1.combine_variants(v2);
	v1.combine_variants(v3);

	REQUIRE(v1.nr_of_alleles() == 3);
	REQUIRE(v1.nr_of_paths() == 4);
	REQUIRE(v1.get_allele_string(0) == "ACTGAGACTG");
	REQUIRE(v1.get_allele_string(1) == "TCTACCACTG");
	REQUIRE(v1.get_allele_string(2) == "TCTACCACTGTC");
	REQUIRE(v1.get_chromosome() == "chr2");
	REQUIRE(v1.get_start_position() == 4);
	REQUIRE(v1.get_end_position() == 14);
	REQUIRE(v1.get_allele_on_path(0) == 0);
	REQUIRE(v1.get_allele_on_path(1) == 0);
	REQUIRE(v1.get_allele_on_path(2) == 2);
	REQUIRE(v1.get_allele_on_path(3) == 1);
}

TEST_CASE("Variant separate_variants", "[Variant separate_variants]") {
	Variant v1("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATT"}, {0,0,1});
	Variant v2("GCT", "CCC", "chr1", 15, 16, {"A", "G"}, {0,1,0});
	Variant v3("ACC", "GGC", "chr1", 18, 19, {"C", "CTA"}, {0,1,1});
	Variant v4("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATT"}, {0,0,1});

	// combine and sepatate 2 variants
	v1.combine_variants(v2);
	vector<Variant> single_variants;
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 2);
	REQUIRE(single_variants[0] == v4);
	REQUIRE(single_variants[1] == v2);

	// combine and separate 3 variants
	v1.combine_variants(v3);
	single_variants.clear();
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 3);
	REQUIRE(single_variants[0] == v4);
	REQUIRE(single_variants[1] == v2);
	REQUIRE(single_variants[2] == v3);

	// add flanking sequences and then separate
	v1.add_flanking_sequence();
	single_variants.clear();
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 3);
	REQUIRE(single_variants[0] == v4);
	REQUIRE(single_variants[1] == v2);
	REQUIRE(single_variants[2] == v3);

	// separate single variant
	single_variants.clear();
	v4.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_variants[0] == v4);
}

TEST_CASE("Variant separate_variants_likelihoods", "Variant separate_variants_likelihoods") {
	Variant v1 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	Variant v2 ("AACT", "ACTG", "chr2", 7, 10, {"GAG", "ACC"}, {0,0,1,1});
	Variant v3 ("GACT", "GGAA", "chr2", 13, 14, {"G", "GTC"}, {0,0,1,0});
	Variant v4 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});

	GenotypingResult g;
	g.add_to_likelihood(0,0,0.05);
	g.add_to_likelihood(0,1,0.05);
	g.add_to_likelihood(1,1,0.0);
	g.add_to_likelihood(0,2,0.3);
	g.add_to_likelihood(1,2,0.5);
	g.add_to_likelihood(2,2,0.1);
	g.add_first_haplotype_allele(0);
	g.add_second_haplotype_allele(2);

	v1.combine_variants(v2);
	v1.combine_variants(v3);
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	v1.separate_variants(&single_variants, &g, &single_genotypes);
	REQUIRE(single_variants.size() == 3);
	REQUIRE(single_variants[0] == v4);
	REQUIRE(single_variants[1] == v2);
	REQUIRE(single_variants[2] == v3);
	REQUIRE(single_genotypes.size() == 3);

	// expected genotype likelihoods
	vector<vector<double>> expected = { {0.05,0.35,0.6}, {0.05,0.35,0.6}, {0.1,0.8,0.1} };
	pair<unsigned char,unsigned char> expected_haplotype = make_pair(0,1);

	// computed genotype likelihoods
	for (size_t i = 0; i < 3; ++i) {
		vector<long double> computed = single_genotypes[i].get_all_likelihoods(2);
		REQUIRE(computed.size() == expected[i].size());
		for (size_t j = 0; j < expected[i].size(); ++j) {
			REQUIRE(doubles_equal(computed[j], expected[i][j]));
		}
		// here, all haplotypes should be 0|1
		REQUIRE(single_genotypes[i].get_haplotype() == expected_haplotype);
	}
}

TEST_CASE("Variant separate_variants_single", "[Variants separate_variants_single]") {
	Variant v ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	GenotypingResult g;
	g.add_to_likelihood(0,0,0.1);
	g.add_to_likelihood(0,1,0.7);
	g.add_to_likelihood(1,1,0.2);

	// separate single variant
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(single_variants[0] == v);

	// same with flanks added
	v.add_flanking_sequence();
	single_variants.clear();
	single_genotypes.clear();
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	v.remove_flanking_sequence();
	REQUIRE(single_variants[0] == v);
}

TEST_CASE("Variant separate_variants_single2", "[Variants separate_variants_single]") {
	Variant v ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {1,1});
	GenotypingResult g;
	g.add_to_likelihood(0,0,0.1);
	g.add_to_likelihood(0,1,0.7);
	g.add_to_likelihood(1,1,0.2);

	// separate single variant
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(single_variants[0] == v);

	// same with flanks added
	v.add_flanking_sequence();
	single_variants.clear();
	single_genotypes.clear();
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	v.remove_flanking_sequence();
	REQUIRE(single_variants[0] == v);
}

TEST_CASE("Variant separate_variants_single3", "[Variants separate_variants_single2]") {
	Variant v ("AAAAAAAAAAAGCCTTTTAACTACTGAAAG", "AAAAAAAAAAAAAAGCACAAGGAAGAAATT", "chr16", 45143, 45144, {"T", "TA"}, {0,0,1,0,0,0,0,0,0,0});
	vector<Variant> single_variants;
	v.add_flanking_sequence();
	v.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 1);
	v.remove_flanking_sequence();
	REQUIRE(single_variants[0] == v);
}

TEST_CASE("Variant uncovered_alleles", "[Variant uncovered_alleles]") {
	Variant v1("AAA", "TCA", "chr1", 4, 5, {"A", "T", "G"}, {0,0});
	Variant v2("AAT", "AAG", "chr1", 6, 7, {"C", "T"}, {0,0});
	Variant v3("CAA", "CCC", "chr1", 9, 10, {"G", "A"}, {0,0});
	Variant v4("AAA", "TCA", "chr1", 4, 5, {"A", "T", "G"}, {0,0});

	REQUIRE(v1.nr_of_alleles() == 3);
	REQUIRE(v2.nr_of_alleles() == 2);
	REQUIRE(v3.nr_of_alleles() == 2);

	// combine variants
	v1.combine_variants(v2);
	v1.combine_variants(v3);
	REQUIRE(v1.nr_of_alleles() == 1);
	REQUIRE(v1.get_allele_string(0) == "ATCAAG");

	// separate variants
	vector<Variant> single_vars;
	v1.separate_variants(&single_vars);

	REQUIRE(single_vars.size() == 3);
	REQUIRE(single_vars[0].nr_of_alleles() == 3);
	REQUIRE(single_vars[1].nr_of_alleles() == 2);
	REQUIRE(single_vars[2].nr_of_alleles() == 2);

	REQUIRE(single_vars[0] == v4);
	REQUIRE(single_vars[1] == v2);
	REQUIRE(single_vars[2] == v3);
}

TEST_CASE("Variant uncovered_single", "[Variant uncovered_single]") {
	Variant v1("AAA", "TTT", "chr1", 5, 6, {"A", "G", "T"}, {0,0,1,0});
	Variant v2("AAA", "TTT", "chr1", 5, 6, {"A", "G", "T"}, {0,0,1,0});
	vector<Variant> single_vars;
	v1.separate_variants(&single_vars);
	REQUIRE(single_vars[0] == v2);
}

TEST_CASE("Variant combine_combined", "[Variant combine_combined]") {
	Variant v1("AAA", "TCA", "chr1", 4, 5, {"A", "T", "G"}, {0,0});
	Variant v2("AAT", "AAG", "chr1", 6, 7, {"C", "T"}, {0,1});
	Variant v3("CAA", "CCC", "chr1", 9, 10, {"G", "A"}, {0,0});
	Variant v4("AAA", "TCA", "chr1", 4, 5, {"A", "T", "G"}, {0,0});
	Variant v5("AAT", "AAG", "chr1", 6, 7, {"C", "T"}, {0,1});

	// combine v2 and v3 first
	v2.combine_variants(v3);
	// combine_variants it with v1
	v1.combine_variants(v2);

	REQUIRE(v1.nr_of_alleles() == 2);
	REQUIRE(v1.get_allele_string(0) == "ATCAAG");
	REQUIRE(v1.get_allele_string(1) == "ATTAAG");

	// separate variants again
	vector<Variant> single_vars;
	v1.separate_variants(&single_vars);
	REQUIRE(single_vars.size() == 3);
	REQUIRE(single_vars[0] == v4);
	REQUIRE(single_vars[1] == v5);
	REQUIRE(single_vars[2] == v3);
}

TEST_CASE("Variant combine_combined2", "[Variant combine_combined]") {
	Variant v1 ("AAA", "TGC", "chr1", 4, 5, {"A", "G"}, {0,0,0,0,0,0,1,0,0,0});
	Variant v2 ("AAT", "CCG", "chr1", 6, 7, {"G", "C"},  {0,0,0,0,0,0,1,0,0,0});
	Variant v3 ("GCC", "GGG", "chr1", 9, 10, {"G", "C"}, {0,0,0,0,0,0,0,1,0,0});
	Variant v4 ("AAA", "TGC", "chr1", 4, 5, {"A", "G"}, {0,0,0,0,0,0,1,0,0,0});

	v1.combine_variants(v2);
	v1.combine_variants(v3);

	REQUIRE(v1.nr_of_alleles() == 3);
	REQUIRE(v1.get_allele_string(0) == "ATGCCG");
	REQUIRE(v1.get_allele_string(1) == "ATGCCC");
	REQUIRE(v1.get_allele_string(2) == "GTCCCG");

	GenotypingResult g;
	g.add_to_likelihood(0,0,0.9);
	g.add_to_likelihood(0,1,0.05);
	g.add_to_likelihood(0,2,0.05);
	g.add_first_haplotype_allele(0);
	g.add_second_haplotype_allele(2);

	vector<Variant> single_vars;
	vector<GenotypingResult> single_genotypes;
	v1.separate_variants(&single_vars, &g, &single_genotypes);
	REQUIRE(single_vars.size() == 3);
	REQUIRE(single_vars[0] == v4);
	REQUIRE(single_vars[1] == v2);
	REQUIRE(single_vars[2] == v3);

	REQUIRE(single_genotypes.size() == 3);
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,0), 0.95));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,1), 0.05));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(1,1), 0.0));

	REQUIRE(single_genotypes[0].get_haplotype() == pair<unsigned char,unsigned char>(0,1));
	REQUIRE(single_genotypes[1].get_haplotype() == pair<unsigned char,unsigned char>(0,1));
	REQUIRE(single_genotypes[2].get_haplotype() == pair<unsigned char,unsigned char>(0,0));
}

TEST_CASE("Variant get_paths_of_allele", "[Variant get_paths_of_allele]") {
	Variant v1("AAA", "TTA", "chr1", 10, 14, {"ATGC", "ATT", "TT"}, {0,1,2});
	vector<size_t> result;
	v1.get_paths_of_allele(0, result);
	REQUIRE(result == vector<size_t>({0}));
	result.clear();
	v1.get_paths_of_allele(1, result);
	REQUIRE(result == vector<size_t>({1}));
	result.clear();
	v1.get_paths_of_allele(2, result);
	REQUIRE(result == vector<size_t>({2}));
	result.clear();

	Variant v2("AAA", "TTA", "chr1", 10, 14, {"ATGC", "ATT"}, {0,1,0,1,1});
	v2.get_paths_of_allele(0, result);
	REQUIRE(result == vector<size_t>({0,2}));
	result.clear();
	v2.get_paths_of_allele(1, result);
	REQUIRE(result == vector<size_t>({1,3,4}));
}
