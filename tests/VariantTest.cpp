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

	vector<unsigned char> path_to_allele = {0, 0, 2, 1};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers> (new UniqueKmers(0, path_to_allele));
	vector<unsigned char> alleles1 = {0};
	vector<unsigned char> alleles2 = {1};
	vector<unsigned char> alleles3 = {2};
	
	for (size_t i = 0; i < 3; ++i) {
		u->insert_kmer(20, alleles1);
	}
	for (size_t i = 0; i < 9; ++i) {
		u->insert_kmer(20, alleles2);
	}
	for (size_t i = 0; i < 2; ++i) {
		u->insert_kmer(20, alleles3);
	}

	v1.combine_variants(v2);
	v1.combine_variants(v3);
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	vector<VariantStats> single_stats;
	v1.separate_variants(&single_variants, &g, &single_genotypes);
	v1.variant_statistics(u, single_stats);
	REQUIRE(single_variants.size() == 3);
	REQUIRE(single_variants[0] == v4);
	REQUIRE(single_variants[1] == v2);
	REQUIRE(single_variants[2] == v3);
	REQUIRE(single_genotypes.size() == 3);
	REQUIRE(single_stats.size() == 3);

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
		REQUIRE(single_stats[i].nr_unique_kmers == 14);
	}

	vector<vector<unsigned int>> expected_counts = { {3,11}, {3,11}, {12,2} };
	vector<vector<string>> expected_alleles = { {"A", "T"}, {"GAG", "ACC"}, {"G", "GTC"} };
	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(single_stats[i].kmer_counts[0] == expected_counts[i][0]);
		REQUIRE(single_stats[i].kmer_counts[1] == expected_counts[i][1]);
		REQUIRE(single_variants[i].get_allele_string(0) == expected_alleles[i][0]);
		REQUIRE(single_variants[i].get_allele_string(1) == expected_alleles[i][1]);
	}
	
}

TEST_CASE("Variant separate_variants_single", "[Variants separate_variants_single]") {
	Variant v ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	GenotypingResult g;
	g.add_to_likelihood(0,0,0.1);
	g.add_to_likelihood(0,1,0.7);
	g.add_to_likelihood(1,1,0.2);

	vector<unsigned char> allele_paths = {0, 0, 1, 1};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers> (new UniqueKmers(0, allele_paths));
	vector<unsigned char> alleles1 = {0,1};
	vector<unsigned char> alleles2 = {1};
	for (size_t i = 0; i < 10; ++i) {
		u->insert_kmer(20, alleles1);
	}
	for (size_t i = 0; i < 20; ++i) {
		u->insert_kmer(30, alleles2);
	}

	// separate single variant
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	vector<VariantStats> single_stats;
	v.separate_variants(&single_variants, &g, &single_genotypes);
	v.variant_statistics(u, single_stats);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,0), 0.1));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,1), 0.7));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(1,1), 0.2));
	REQUIRE(single_variants[0] == v);
	REQUIRE(single_stats[0].kmer_counts[0] == 10);
	REQUIRE(single_stats[0].kmer_counts[1] == 30);

	// same with flanks added
	v.add_flanking_sequence();
	single_variants.clear();
	single_genotypes.clear();
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,0), 0.1));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,1), 0.7));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(1,1), 0.2));
	v.remove_flanking_sequence();
	REQUIRE(single_variants[0] == v);
	REQUIRE(single_stats[0].nr_unique_kmers == 30);
	REQUIRE(single_stats[0].kmer_counts[0] == 10);
	REQUIRE(single_stats[0].kmer_counts[1] == 30);
}

TEST_CASE("Variant separate_variants_single2", "[Variants separate_variants_single]") {
	Variant v ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {1,1});
	GenotypingResult g;
	g.add_to_likelihood(0,0,0.1);
	g.add_to_likelihood(0,1,0.7);
	g.add_to_likelihood(1,1,0.2);

	vector<unsigned char> path_to_allele = {1, 1};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers> (new UniqueKmers (0, path_to_allele));
	vector<unsigned char> alleles1 = {0};
	vector<unsigned char> alleles2 = {1};
	u->insert_kmer(20, alleles1);
	u->insert_kmer(20, alleles1);
	u->insert_kmer(20, alleles2);
	u->insert_kmer(20, alleles2);
	u->insert_kmer(20, alleles2);
	u->insert_kmer(20, alleles2);

	// separate single variant
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	vector<VariantStats> single_stats;
	v.separate_variants(&single_variants, &g, &single_genotypes);
	v.variant_statistics(u, single_stats);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(single_variants[0] == v);
	REQUIRE(single_stats[0].nr_unique_kmers == 6);
	// allele 0 is not covered thus the kmer count should be -1
	REQUIRE(single_stats[0].kmer_counts[0] == -1);
	REQUIRE(single_stats[0].kmer_counts[1] == 4);

	// same with flanks added
	v.add_flanking_sequence();
	single_variants.clear();
	single_genotypes.clear();
	v.separate_variants(&single_variants, &g, &single_genotypes);

	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	v.remove_flanking_sequence();
	REQUIRE(single_variants[0] == v);
	REQUIRE(single_stats[0].nr_unique_kmers == 6);
	// allele 0 is not covered thus the kmer count should be -1
	REQUIRE(single_stats[0].kmer_counts[0] == -1);
	REQUIRE(single_stats[0].kmer_counts[1] == 4);
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

TEST_CASE("Variant nr_missing_alleles", "[Variant nr_missing_alleles]") {
	Variant v1("AAA", "TTT", "chr1", 5, 6, {"A", "GNN", "T"}, {0,1,1,2});
	REQUIRE(v1.nr_missing_alleles() == 2);

	Variant v2("AAAN", "TTTN", "chr1", 5, 6, {"A", "G", "T"}, {0,0,1,0});
	REQUIRE(v2.nr_missing_alleles() == 0);
	v2.add_flanking_sequence();
	REQUIRE(v2.nr_missing_alleles() == 4);
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
	Variant v1 ("AAA", "TGC", "chr1", 4, 5, {"A", "G"},  {0,0,0,0,0,0,1,0,0,0});
	Variant v2 ("AAT", "CCG", "chr1", 6, 7, {"G", "C"},  {0,0,0,0,0,0,1,0,0,0});
	Variant v3 ("GCC", "GGG", "chr1", 9, 10, {"G", "C"}, {0,0,0,0,0,0,0,1,0,0});
	Variant v4 ("AAA", "TGC", "chr1", 4, 5, {"A", "G"},  {0,0,0,0,0,0,1,0,0,0});

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

	vector<unsigned char> path_to_allele = {0, 0, 0, 0, 0, 0, 2, 1, 0};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers>( new UniqueKmers(0, path_to_allele));
	vector<unsigned char> alleles1 = {0};
	vector<unsigned char> alleles2 = {1};
	vector<unsigned char> alleles3 = {2};	
	for (size_t i = 0; i < 10; i++) {
		u->insert_kmer(20, alleles1);
	}
	for (size_t i = 0; i < 2; ++i) {
		u->insert_kmer(30, alleles2);
	}
	for (size_t i = 0; i < 4; ++i) {
		u->insert_kmer(25, alleles3);
	}

	vector<Variant> single_vars;
	vector<GenotypingResult> single_genotypes;
	vector<VariantStats> single_stats;
	v1.separate_variants(&single_vars, &g, &single_genotypes);
	v1.variant_statistics(u, single_stats);
	REQUIRE(single_vars.size() == 3);
	REQUIRE(single_vars[0] == v4);
	REQUIRE(single_vars[1] == v2);
	REQUIRE(single_vars[2] == v3);
	REQUIRE(single_stats.size() == 3);

	REQUIRE(single_genotypes.size() == 3);
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,0), 0.95));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(0,1), 0.05));
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(1,1), 0.0));

	REQUIRE(single_genotypes[0].get_haplotype() == pair<unsigned char,unsigned char>(0,1));
	REQUIRE(single_genotypes[1].get_haplotype() == pair<unsigned char,unsigned char>(0,1));
	REQUIRE(single_genotypes[2].get_haplotype() == pair<unsigned char,unsigned char>(0,0));

	vector<vector<unsigned int>> expected_counts = { {12,4}, {12,4}, {14,2} };
	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(single_stats[i].kmer_counts[0] == expected_counts[i][0]);
		REQUIRE(single_stats[i].kmer_counts[1] == expected_counts[i][1]);
	}
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

TEST_CASE("Variant allele_frequency", "Variant allele_frequency") {
	Variant v1("AAA", "TTA", "chr1", 10, 14, {"ATGC", "ATT", "TT"}, {0,1,2});
	REQUIRE ( doubles_equal(v1.allele_frequency(0), 1.0/3.0) );
	REQUIRE ( doubles_equal(v1.allele_frequency(1), 1.0/3.0) );
	REQUIRE ( doubles_equal(v1.allele_frequency(2), 1.0/3.0) );
	REQUIRE ( doubles_equal(v1.allele_frequency(0, true), 0.0) );
	REQUIRE ( doubles_equal(v1.allele_frequency(1, true), 0.5) );
	REQUIRE ( doubles_equal(v1.allele_frequency(2, true), 0.5) );

	Variant v2("AAA", "TGC", "chr1", 4, 5, {"A", "G"}, {0,0,0,0,0,0,1,0,0,0});
	REQUIRE ( doubles_equal(v2.allele_frequency(0), 0.9) );
	REQUIRE ( doubles_equal(v2.allele_frequency(1), 0.1) );
	REQUIRE ( doubles_equal(v2.allele_frequency(0, true), 8.0/9.0) );
	REQUIRE ( doubles_equal(v2.allele_frequency(1, true), 1.0/9.0) );

	Variant v3("AAA", "TGC", "chr1", 4, 5, {"A", "G", "C"}, {0,0,1,0,2,0,1,0,0,0});
	REQUIRE ( doubles_equal(v3.allele_frequency(0), 0.7) );
	REQUIRE ( doubles_equal(v3.allele_frequency(1), 0.2) );
	REQUIRE ( doubles_equal(v3.allele_frequency(2), 0.1) );
	REQUIRE ( doubles_equal(v3.allele_frequency(0, true), 6.0/9.0) );
	REQUIRE ( doubles_equal(v3.allele_frequency(1, true), 2.0/9.0) );
	REQUIRE ( doubles_equal(v3.allele_frequency(2, true), 1.0/9.0) );
}

TEST_CASE("Variant separate_variants_likelihoods_uncovered", "Variant separate_variants_likelihoods_uncovered") {
	Variant v1 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,1}); //, "VAR1");
	// second allele is not covered by any path
	Variant v2 ("AACT", "ACTG", "chr2", 7, 8, {"G", "C", "T"}, {0,2}); //, "VAR2");

	GenotypingResult g;
	g.add_to_likelihood(0,0,0.05);
	g.add_to_likelihood(0,1,0.05);
	g.add_to_likelihood(1,1,0.9);
	g.add_first_haplotype_allele(0);
	g.add_second_haplotype_allele(0);

	vector<unsigned char> path_to_allele = {0, 1};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers> (new UniqueKmers (0, path_to_allele));
	vector<unsigned char> alleles1 = {0};
	vector<unsigned char> alleles2 = {1};
	for (size_t i = 0; i < 3; ++i) {
		u->insert_kmer(20, alleles1);
	}
	for (size_t i = 0; i < 9; ++i) {
		u->insert_kmer(20, alleles2);
	}

	v1.combine_variants(v2);
	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	vector<VariantStats> variant_stats;
	REQUIRE(v1.get_id() == ".");
	v1.separate_variants(&single_variants, &g, &single_genotypes);
	v1.variant_statistics(u, variant_stats);
	REQUIRE(single_variants.size() == 2);
	REQUIRE(single_genotypes.size() == 2);
	REQUIRE(variant_stats.size() == 2);

	// expected genotype likelihoods
	// order of alleles changes: uncovered allele has index 2 after separation
	vector<vector<double>> expected = { {0.05,0.05,0.9}, {0.05, 0.0, 0.0, 0.05, 0.0, 0.9}};
	pair<unsigned char,unsigned char> expected_haplotype = make_pair(0,0);
	vector<unsigned int> nr_alleles = {2,3};

	// computed genotype likelihoods
	for (size_t i = 0; i < 2; ++i) {
		vector<long double> computed = single_genotypes[i].get_all_likelihoods(nr_alleles[i]);
		REQUIRE(computed.size() == expected[i].size());
		for (size_t j = 0; j < expected[i].size(); ++j) {
			REQUIRE(doubles_equal(computed[j], expected[i][j]));
		}
		// here, all haplotypes should be 0|0
		REQUIRE(single_genotypes[i].get_haplotype() == expected_haplotype);
		REQUIRE(variant_stats[i].nr_unique_kmers == 12);
	}

	// uncovered allele should have count -1
	vector<vector<int>> expected_counts = { {3,9}, {3,-1,9} };
	vector<vector<string>> expected_alleles = { {"A", "T"}, {"G", "C", "T"} };

	REQUIRE(variant_stats[0].kmer_counts[0] == expected_counts[0][0]);
	REQUIRE(variant_stats[0].kmer_counts[1] == expected_counts[0][1]);

	REQUIRE(single_variants[0].get_allele_string(0) == expected_alleles[0][0]);
	REQUIRE(single_variants[0].get_allele_string(1) == expected_alleles[0][1]);

	REQUIRE(variant_stats[1].kmer_counts[0] == expected_counts[1][0]);
	REQUIRE(variant_stats[1].kmer_counts[1] == expected_counts[1][1]);
	REQUIRE(variant_stats[1].kmer_counts[2] == expected_counts[1][2]);
	REQUIRE(single_variants[1].get_allele_string(0) == expected_alleles[1][0]);
	REQUIRE(single_variants[1].get_allele_string(1) == expected_alleles[1][1]);
	REQUIRE(single_variants[1].get_allele_string(2) == expected_alleles[1][2]);
}

TEST_CASE("Variant separate_variants_likelihoods_single_uncovered", "[Variant separate_variants_likelihoods_single_uncovered]") {
	Variant v ("ATGA", "CTGA", "chr1", 7, 8, {"A", "T"}, {1,1});
	GenotypingResult g;
	g.add_to_likelihood(1,1,1.0);
	g.add_first_haplotype_allele(1);
	g.add_second_haplotype_allele(1);
	
	vector<unsigned char> path_to_allele = {1, 1};
	shared_ptr<UniqueKmers> u = shared_ptr<UniqueKmers> (new UniqueKmers(0, path_to_allele));
	vector<unsigned char> alleles1 = {0};
	vector<unsigned char> alleles2 = {1};
	u->insert_kmer(20, alleles1);
	u->insert_kmer(30, alleles2);
	u->insert_kmer(25, alleles2);
	u->insert_kmer(20, alleles2);

	vector<Variant> single_variants;
	vector<GenotypingResult> single_genotypes;
	v.separate_variants(&single_variants, &g, &single_genotypes);
	REQUIRE(single_variants.size() == 1);
	REQUIRE(single_genotypes.size() == 1);
	REQUIRE(doubles_equal(single_genotypes[0].get_genotype_likelihood(1,1), 1.0));
	REQUIRE(single_variants[0].get_allele_string(0) == "A");
	REQUIRE(single_variants[0].get_allele_string(1) == "T");

	vector<VariantStats> stats;
	v.variant_statistics(u, stats);
	REQUIRE(stats.size() == 1);
	REQUIRE(stats[0].nr_unique_kmers == 4);
	// allele 0 is not covered by any path and its kmer count should therefore be -1 
	REQUIRE(stats[0].kmer_counts[0] == -1);
	REQUIRE(stats[0].kmer_counts[1] == 3);
	REQUIRE(stats[0].coverage == 0);
}

TEST_CASE("Variant get_id", "[Variant get_id]") {
	Variant v1 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,1}); //, "VAR1");
	Variant v2 ("AACT", "ACTG", "chr2", 7, 8, {"G", "C", "T"}, {0,2}); //, "VAR2");

	REQUIRE(v1.get_id() == ".");
	REQUIRE(v2.get_id() == ".");
	v1.combine_variants(v2);

	REQUIRE(v1.get_id() == ".");
	vector<Variant> single_variants;
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 2);
	REQUIRE(single_variants[0].get_id() == ".");
	REQUIRE(single_variants[1].get_id() == ".");
}

TEST_CASE("Variant get_id2", "[Variant get_id2]") {
	Variant v1 ("AAA", "TGC", "chr1", 4, 5, {"A", "G"}, {0,0,0,0,0,0,1,0,0,0}); //, "VAR1");
	Variant v2 ("AAT", "CCG", "chr1", 6, 7, {"G", "C"},  {0,0,0,0,0,0,1,0,0,0}); //, ".");
	Variant v3 ("GCC", "GGG", "chr1", 9, 10, {"G", "C"}, {0,0,0,0,0,0,0,1,0,0}); //, "VAR2;VAR3");

	REQUIRE(v1.get_id() == ".");
	REQUIRE(v2.get_id() == ".");
	REQUIRE(v3.get_id() == ".");
	v1.combine_variants(v2);
	v1.combine_variants(v3);

	REQUIRE(v1.get_id() == ".");
	vector<Variant> single_variants;
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 3);
	REQUIRE(single_variants[0].get_id() == ".");
	REQUIRE(single_variants[1].get_id() == ".");
	REQUIRE(single_variants[2].get_id() == ".");
}

TEST_CASE("Variant is_undefined_allele", "[Variant is_undefined_allele]"){
	Variant v1("AAN", "TAC", "chr1", 10, 14, {"ATGC", "ATT"}, {0,1});

	REQUIRE(!v1.is_undefined_allele(0));
	REQUIRE(!v1.is_undefined_allele(1));

	v1.add_flanking_sequence();
	REQUIRE(!v1.is_undefined_allele(0));
	REQUIRE(!v1.is_undefined_allele(1));

	Variant v2("GCT", "CCC", "chr1", 15, 17, {"AN", "G"}, {1,0});
	REQUIRE(v2.is_undefined_allele(0));
	REQUIRE(!v2.is_undefined_allele(1));

	v2.add_flanking_sequence();
	REQUIRE(v2.is_undefined_allele(0));
	REQUIRE(!v2.is_undefined_allele(1));
}

TEST_CASE("Variant combine_variants_undefined_flanks", "[Variant combine_variants_undefined_flanks]") {
	Variant v1 ("ATGA", "CNGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	Variant v2 ("AACN", "ACTG", "chr2", 7, 10, {"GAG", "ACC"}, {0,0,1,1});
	Variant v3 ("GACT", "GGAA", "chr2", 13, 14, {"G", "GTC"}, {0,0,1,0});

	v1.combine_variants(v2);
	v1.combine_variants(v3);

	REQUIRE(v1.nr_of_alleles() == 3);
	REQUIRE(v1.nr_of_paths() == 4);
	REQUIRE(v1.get_allele_string(0) == "ACNGAGACTG");
	REQUIRE(v1.get_allele_string(1) == "TCNACCACTG");
	REQUIRE(v1.get_allele_string(2) == "TCNACCACTGTC");
	REQUIRE(v1.get_chromosome() == "chr2");
	REQUIRE(v1.get_start_position() == 4);
	REQUIRE(v1.get_end_position() == 14);
	REQUIRE(v1.get_allele_on_path(0) == 0);
	REQUIRE(v1.get_allele_on_path(1) == 0);
	REQUIRE(v1.get_allele_on_path(2) == 2);
	REQUIRE(v1.get_allele_on_path(3) == 1);

	v1.add_flanking_sequence();
	REQUIRE(!v1.is_undefined_allele(0));
	REQUIRE(!v1.is_undefined_allele(1));
	REQUIRE(!v1.is_undefined_allele(2));

	v1.remove_flanking_sequence();
	REQUIRE(!v1.is_undefined_allele(0));
	REQUIRE(!v1.is_undefined_allele(1));
	REQUIRE(!v1.is_undefined_allele(2));

	v1.add_flanking_sequence();
	vector<Variant> single_variants;
	v1.separate_variants(&single_variants);

	for (auto v : single_variants) {
		REQUIRE(!v.is_undefined_allele(0));
		REQUIRE(!v.is_undefined_allele(1));
	}
}

TEST_CASE("Variant combine_variants_undefined_alleles", "[Variant combine_variants_undefined_flanks]") {
	Variant v1 ("ATGA", "CTGA", "chr2", 4, 5, {"A", "T"}, {0,0,1,1});
	Variant v2 ("AACT", "ACTG", "chr2", 7, 10, {"GNG", "ACC"}, {0,0,1,1});

	v1.combine_variants(v2);
	REQUIRE(v1.nr_of_alleles() == 2);
	REQUIRE(v1.get_allele_string(0) == "ACTGNG");
	REQUIRE(v1.get_allele_string(1) == "TCTACC");

	REQUIRE(v1.is_undefined_allele(0));
	REQUIRE(!v1.is_undefined_allele(1));

	v1.add_flanking_sequence();
	vector<Variant> single_variants;
	v1.separate_variants(&single_variants);
	REQUIRE(single_variants.size() == 2);
	REQUIRE(!single_variants[0].is_undefined_allele(0));
	REQUIRE(!single_variants[0].is_undefined_allele(1));
	REQUIRE(single_variants[1].is_undefined_allele(0));
	REQUIRE(!single_variants[1].is_undefined_allele(1));
}

TEST_CASE("Variant separate_variants_identical", "[Variant separate_variants_identical]") {

	Variant v1("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATGC"}, {0,0,1});
	Variant v2("GCT", "CCN", "chr1", 15, 16, {"A", "A"}, {0,1,0});
	Variant v3("ACC", "GGC", "chr1", 18, 19, {"N", "N"}, {0,1,1});
	Variant v4("AAA", "TAC", "chr1", 10, 14, {"ATGC", "ATGC"}, {0,0,1});


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