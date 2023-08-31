#include "catch.hpp"
#define private public
#include "../src/graphbuilder.hpp"
#include "../src/graph.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/variant.hpp"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <algorithm> 
#include <random>


using namespace std;

TEST_CASE("GraphBuilder get_allele_string", "[GraphBuilder get_allele_string]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/small1-segments.fa", 10, true);

	vector<string> chromosomes;
	v.get_chromosomes(&chromosomes);
	REQUIRE(chromosomes.size() == 2);
	REQUIRE(chromosomes[0] == "chrA");
	REQUIRE(chromosomes[1] == "chrB");
	REQUIRE(v.get_kmer_size() == 10);
	// there should be 4+1 paths (since the reference path is added)
	REQUIRE(v.nr_of_paths() == 5);

	REQUIRE(graph.size() == 2);

	// double check if kmer sizes and chromosomes have been set correctly
	REQUIRE(graph.at("chrA")->get_kmer_size() == 10);
	REQUIRE(graph.at("chrB")->get_kmer_size() == 10);
	REQUIRE(graph.at("chrA")->get_chromosome() == "chrA");
	REQUIRE(graph.at("chrB")->get_chromosome() == "chrB");

	// check if graph segments have been properly constructed
	REQUIRE(graph.at("chrA")->size() == 7);
	REQUIRE(graph.at("chrB")->size() == 2);
	REQUIRE(graph.at("chrA")->get_variant(2).nr_of_alleles() == 3);

	REQUIRE(graph.at("chrA")->get_variant(0).get_allele_string(0) == "GGAATTCCGACATAAGTTA");
	REQUIRE(graph.at("chrA")->get_variant(0).get_allele_string(1) == "GGAATTCCGTCATAAGTTA");

	REQUIRE(graph.at("chrA")->get_variant(1).get_allele_string(0) == "CCTTAGCTACGAAGCCAGT");
	REQUIRE(graph.at("chrA")->get_variant(1).get_allele_string(1) == "CCTTAGCTAGGGGGAAGCCAGT");

	REQUIRE(graph.at("chrA")->get_variant(2).get_allele_string(0) == "GAAGCCAGTGCCCCGAGACGGCCAAA");
	REQUIRE(graph.at("chrA")->get_variant(2).get_allele_string(1) == "GAAGCCAGTTCCCCGAGACGGCCAAA");
	REQUIRE(graph.at("chrA")->get_variant(2).get_allele_string(2) == "GAAGCCAGTTCCCCTACGGCCAAA");
	REQUIRE(graph.at("chrA")->get_variant(2).nr_of_paths() == 5);

	REQUIRE(graph.at("chrA")->get_variant(3).get_allele_string(0) == "ACGTCCGTTCAGCCTTAGC");
	REQUIRE(graph.at("chrA")->get_variant(3).get_allele_string(1) == "ACGTCCGTTTAGCCTTAGC");

	REQUIRE(graph.at("chrA")->get_variant(4).get_allele_string(0) == "CCGATTTTCTTGTGCTATA");
	REQUIRE(graph.at("chrA")->get_variant(4).get_allele_string(1) == "CCGATTTTCCTGTGCTATA");

	REQUIRE(graph.at("chrA")->get_variant(5).get_allele_string(0) == "GGAGGGTATGAAGCCATCAC");
	REQUIRE(graph.at("chrA")->get_variant(5).get_allele_string(1) == "GGAGGGTATTCAGCCATCAC");

	REQUIRE(graph.at("chrA")->get_variant(6).get_allele_string(0) == "TGTGGACTTATTTGGCTAA");
	REQUIRE(graph.at("chrA")->get_variant(6).get_allele_string(1) == "TGTGGACTTGTTTGGCTAA");

	REQUIRE(graph.at("chrB")->get_variant(0).get_allele_string(0) == "CCACTTCATCAAGACACAA");
	REQUIRE(graph.at("chrB")->get_variant(1).get_allele_string(0) == "GAGTATTTTGATCATAAAT");
}



TEST_CASE("GraphBuilder get_overhang", "[GraphBuilder get_overhang]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/small1-segments.fa", 10, false);

	DnaSequence left_overhang;
	DnaSequence right_overhang;
	REQUIRE(graph.size() == 2);
	graph.at("chrA")->get_left_overhang(0, 20, left_overhang);
	graph.at("chrA")->get_right_overhang(0, 20, right_overhang);

	REQUIRE(left_overhang.to_string() == "TTTGGTGATCTGGAATTCCG");
	REQUIRE(right_overhang.to_string() == "CATAAGTTATGCTAAAAAAT");

	graph.at("chrA")->get_left_overhang(1, 20, left_overhang);
	graph.at("chrA")->get_right_overhang(1, 20, right_overhang);
	REQUIRE(left_overhang.to_string() == "GTCTGTTAAGACCTTAGCTA");
	REQUIRE(right_overhang.to_string() == "GAAGCCAGT");

	graph.at("chrA")->get_left_overhang(2, 20, left_overhang);
	graph.at("chrA")->get_right_overhang(2, 20, right_overhang);
	REQUIRE(left_overhang.to_string() == "GAAGCCAGT");
	REQUIRE(right_overhang.to_string() == "ACGGCCAAAACATACCATTT");
	REQUIRE(v.nr_of_paths() == 4);
}

TEST_CASE("GraphBuilder write_path_segments", "[GrpahBuilder write_path_segments]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF file
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/small1-segments.fa", 10, false);

	// compare reference segments to expected sequences
	vector<string> computed = {};
	vector<string> expected = {};

	// read expected reference segments from file
	ifstream expected_sequences("../tests/data/small1-expected-ref-segments.fa");
	string line;
	while (getline(expected_sequences, line)) {
		size_t start = line.find_first_not_of(" \t\r\n");
		size_t end = line.find_last_not_of(" \t\r\n");
		line = line.substr(start, end-start+1);
		expected.push_back(line);
	}

	// read computed reference segments from file
	bool read_next = false;
	ifstream computed_sequences("../tests/data/small1-segments.fa");
	while (getline(computed_sequences, line)) {
		if (line.size() == 0) continue;
		if (line[0] == '>') {
			if (line.find("reference") != string::npos) {
				read_next = true;
			} else {
				read_next = false;
			}
			continue;
		}
		size_t start = line.find_first_not_of(" \t\r\n");
		size_t end = line.find_last_not_of(" \t\r\n");
		line = line.substr(start, end-start+1);
		if (read_next) computed.push_back(line);
	}

	REQUIRE(expected.size() == computed.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(expected[i] == computed[i]);
	}
}

TEST_CASE("GraphBuilder write_path_segments_no_variants", "[GraphBuilder write_path_segments]") {
	string vcf = "../tests/data/empty.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, false);

	vector<string> computed = {};
	vector<string> expected = {};

	// since no variants are given, expected sequences are all those in fasta
	string line;
	ifstream expected_sequences(fasta);
	string sequence = "";
	while (getline(expected_sequences, line)) {
		if (line.size() == 0) continue;
		if (line[0] == '>') {
			if (sequence != "") {
				transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
				expected.push_back(sequence);
				sequence = "";
			}
			continue;
		}
		size_t start = line.find_first_not_of(" \t\r\n");
		size_t end = line.find_last_not_of(" \t\r\n");
		line = line.substr(start, end-start+1);
		sequence += line;
	}
	expected.push_back(sequence);

	// read computed reference segments from file
	bool read_next = false;
	ifstream computed_sequences("../tests/data/empty-segments.fa");
	while (getline(computed_sequences, line)) {
		if (line.size() == 0) continue;
		if (line[0] == '>') {
			continue;
		}
		size_t start = line.find_first_not_of(" \t\r\n");
		size_t end = line.find_last_not_of(" \t\r\n");
		line = line.substr(start, end-start+1);
		computed.push_back(line);
	}

	REQUIRE(expected.size() == computed.size());
	for (size_t i = 0; i < expected.size(); ++i) {
		REQUIRE(expected[i] == computed[i]);
	}
}

TEST_CASE("GraphBuilder write_genotypes_of", "[GraphBuilder write_genotypes_of]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF file
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, false);

	vector<string> chromosomes;
	vector<string> expected_chromosomes = {"chrA", "chrB"};
	v.get_chromosomes(&chromosomes);

	REQUIRE(chromosomes.size() == expected_chromosomes.size());
	REQUIRE(chromosomes[0] == expected_chromosomes[0]);
	REQUIRE(chromosomes[1] == expected_chromosomes[1]);

	// generate a GenotypingResult for chrA
	vector<GenotypingResult> genotypes_chrA(7);
	
	for (size_t i = 0; i < 7; ++i) {
		if (i == 2) continue;
		GenotypingResult r;
		r.add_to_likelihood(0,0,0.2);
		r.add_to_likelihood(0,1,0.7);
		r.add_to_likelihood(1,1,0.1);
		genotypes_chrA[i] = r;
	}

	// third variant is multiallelic
	GenotypingResult r;
	r.add_to_likelihood(0,0,0.2);
	r.add_to_likelihood(0,1,0.0);
	r.add_to_likelihood(0,2,0.2);
	r.add_to_likelihood(1,1,0.0);
	r.add_to_likelihood(1,2,0.5);
	r.add_to_likelihood(2,2,0.1);
	r.add_first_haplotype_allele(2);
	r.add_second_haplotype_allele(1);
	genotypes_chrA[2] = r;


	// generate a GenotypingResult for chrB
	vector<GenotypingResult> genotypes_chrB(2);
	for (size_t i = 0; i < 2; ++i) {
		GenotypingResult r;
		r.add_to_likelihood(0,0,0.1);
		r.add_to_likelihood(0,1,0.1);
		r.add_to_likelihood(1,1,0.8);
		genotypes_chrB[i] = r;
	}

	graph.at("chrA")->write_genotypes("../tests/data/small1-genotypes-graph.vcf", genotypes_chrA, true, "HG0");
	graph.at("chrB")->write_genotypes("../tests/data/small1-genotypes-graph.vcf", genotypes_chrB, false, "HG0");

	graph.at("chrA")->write_phasing("../tests/data/small1-phasing-graph.vcf", genotypes_chrA, true, "HG0");
	graph.at("chrB")->write_phasing("../tests/data/small1-phasing-graph.vcf", genotypes_chrB, false, "HG0");
}

TEST_CASE("GraphBuilder broken_vcfs", "[GraphBuilder broken_vcfs]") {
	string no_paths = "../tests/data/no-paths.vcf";
	string malformatted = "../tests/data/malformatted-vcf1.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;

	CHECK_THROWS(GraphBuilder(no_paths, fasta, graph, "../tests/data/empty-segments.fa", 10, false));
	CHECK_THROWS(GraphBuilder(malformatted, fasta, graph, "../tests/data/empty-segments.fa", 10, false));
}

TEST_CASE("GraphBuilder no-alt-alleles", "[GraphBuilder no-alt-alleles]") {
	string vcf = "../tests/data/no-alt-alleles.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	GraphBuilder(vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, false);

	// should have skipped variant for which no alt alleles are given
	REQUIRE(graph.at("chrA")->size() == 1);
}

TEST_CASE("GraphBuilder overlapping variants", "[GraphBuilder overlapping variants]") {
	string vcf = "../tests/data/overlapping-variants.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	REQUIRE_THROWS(GraphBuilder(vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, false));
}

TEST_CASE("GraphBuilder get_chromosomes", "[GraphBuilder get_chromosomes]") {
	string vcf1 = "../tests/data/small1.vcf";
	string vcf2 = "../tests/data/small2.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v1(vcf1, fasta, graph, "../tests/data/empty-segments.fa", 10, false);

	vector<string> chromosomes;
	v1.get_chromosomes(&chromosomes);
	vector<string> expected1 = {"chrA", "chrB"};
	REQUIRE(chromosomes == expected1);

	graph.clear();
	GraphBuilder v2(vcf2, fasta, graph, "../tests/data/empty-segments.fa", 10, false);
	chromosomes.clear();
	v2.get_chromosomes(&chromosomes);
	vector<string> expected2 = {"chrB", "chrC", "chrA"};
	REQUIRE(chromosomes == expected2);
}

TEST_CASE("GraphBuilder construct_index", "[GraphBuilder construct_index]") {
	vector<string> sequences = {"TTTTT", "AATAGTAAAGTTATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	vector<DnaSequence> alleles;
	for (auto s : sequences) {
		alleles.push_back(DnaSequence(s));
	}
	vector<unsigned char> expected = {1,0,2,3};
	REQUIRE(graph_construct_index(alleles, true) == expected);
}

TEST_CASE("GraphBuilder variant_ids1", "[GraphBuilder variant_ids1]") {
	string vcf = "../tests/data/small1-ids.vcf";
	string fasta = "../tests/data/small1.fa";

	Graph graph;

	vector<string> sequences_ref = {"TGGG", "AATAGTAAAGTTATA", "GTAGATAGATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	vector<string> sequences = {"AATAGTAAAGTTATA", "GTAGATAGATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	map<string,string> sequence_to_id = {{"AATAGTAAAGTTATA", "var1"}, {"GTAGATAGATA", "var2"}, {"AATAGTAAAGTGATA", "var3"}, {"GGGTG", "var4"}, {"TTG", "var5:var6"}};

	vector<DnaSequence> alleles;
	for (auto s : sequences_ref) {
		alleles.push_back(DnaSequence(s));
	}
	string chromosome = "chrA";
	vector<string> variant_ids = {"var1", "var2", "var3", "var4", "var5:var6"};

	graph.insert_ids(alleles, variant_ids, true);

	for (size_t i = 0; i < 10; ++i) {
		// shuffle alleles
		auto rng = std::default_random_engine {};
		std::shuffle(std::begin(sequences), std::end(sequences), rng);
		string expected = "";
		size_t index = 0;
		for (auto s : sequences) {
			if (index > 0) expected += ',';
			expected += sequence_to_id[s];
			index += 1;
		}
		string result = graph.get_ids(sequences, 0 , false);
		REQUIRE(result == expected);
	}
}

TEST_CASE("GraphBuilder variant_ids2", "[GraphBuilder variant_ids2]") {
	string vcf = "../tests/data/small1-ids.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, true);
	vector<GenotypingResult> genotypes(2);
	REQUIRE(graph.size() == 1);
	graph.at("chrA")->write_genotypes("../tests/data/small1-ids-genotypes.vcf", genotypes, true, "sample");
}

TEST_CASE("GraphBuilder close_to_start", "[GraphBuilder close_to_start]") {
	string vcf = "../tests/data/close.vcf";
	string fasta = "../tests/data/close.fa";
	vector<GenotypingResult> genotypes(1);

	map<string, shared_ptr<Graph>> graph;
	GraphBuilder v(vcf, fasta, graph, "../tests/data/empty-segments.fa", 31, true);
	graph.at("chr10")->write_genotypes("../tests/data/small1-ids-close-graph.vcf", genotypes, true, "sample");
}

TEST_CASE("GraphBuilder non_existing_path", "[GraphBuilder non_existing_path]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	map<string, shared_ptr<Graph>> graph;
	CHECK_THROWS(GraphBuilder (vcf, fasta, graph, "nonexistent/paths_segments.fasta", 10, false));
}

TEST_CASE("GraphBuilder too_large_panel", "[GraphBuilder too_large_panel]") {
	string vcf = "../tests/data/large-panel.vcf";
	string fasta = "../tests/data/small1.fa";
	// there are more than 256 paths in the VCF, the implementation cannot handle this and should throw an error
	map<string, shared_ptr<Graph>> graph;
	CHECK_THROWS(GraphBuilder (vcf, fasta, graph, "nonexistent/paths_segments.fasta", 10, false));
	CHECK_THROWS(GraphBuilder (vcf, fasta, graph, "nonexistent/paths_segments.fasta", 10, true));
}

TEST_CASE("GraphBuilder too_many_alleles", "[GraphBuilder too_many_alleles]") {
	string vcf = "../tests/data/many-alleles.vcf";
	string fasta = "../tests/data/small1.fa";
	// there are more than 256 alleles in the VCF, the implementation cannot handle this and should throw an error
	map<string, shared_ptr<Graph>> graph;
	CHECK_THROWS(GraphBuilder (vcf, fasta, graph, "nonexistent/paths_segments.fasta", 10, false));
}

TEST_CASE("GraphBuilder unknown_alleles", "[GraphBuilder unknown_alleles]") {
	string vcf = "../tests/data/small3.vcf";
	string fasta = "../tests/data/small1.fa";
	// vcf contains unknown alleles in panel (".")
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder (vcf, fasta, graph, "../tests/data/empty-segments.fa", 10, false);	
}

TEST_CASE("GraphBuilder unknown_alleles2", "[GraphBuilder unknown_alleles2]") {
	vector<string> alleles = {"G", "A", "C", "AAA"};
	shared_ptr<Variant> v1 (new Variant ("AAAA", "TTTT", "chr1", 10, 11, {"G", "AAA", "CN", "C", "N", "A"}, {0,1,2}));
	Graph g;
	vector<shared_ptr<Variant>> cluster = {v1};
	vector<vector<string>> variant_ids = { {"var1", "var2", "var3"} };
	g.add_variant_cluster(&cluster, variant_ids, true);

	string computed_ids = g.get_ids(alleles, 0, true);
	string expected_ids = "var3,var2,var1";
	REQUIRE(expected_ids == computed_ids);
}