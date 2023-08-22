#include "catch.hpp"
#define private public
#include "../src/variantreader.hpp"
#include "../src/uniquekmers.hpp"
#include <vector>
#include <string>
#include <algorithm> 
#include <random>


using namespace std;

TEST_CASE("VariantReader get_allele_string", "[VariantReader get_allele_string]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	VariantReader v(vcf, fasta, 10, true);
	REQUIRE(v.get_kmer_size() == 10);
	REQUIRE(v.size_of("chrA") == 7);
	REQUIRE(v.size_of("chrB") == 2);
	REQUIRE(v.get_variant("chrA", 2).nr_of_alleles() == 3);

	REQUIRE(v.get_variant("chrA", 0).get_allele_string(0) == "GGAATTCCGACATAAGTTA");
	REQUIRE(v.get_variant("chrA", 0).get_allele_string(1) == "GGAATTCCGTCATAAGTTA");

	REQUIRE(v.get_variant("chrA", 1).get_allele_string(0) == "CCTTAGCTACGAAGCCAGT");
	REQUIRE(v.get_variant("chrA", 1).get_allele_string(1) == "CCTTAGCTAGGGGGAAGCCAGT");

	REQUIRE(v.get_variant("chrA", 2).get_allele_string(0) == "GAAGCCAGTGCCCCGAGACGGCCAAA");
	REQUIRE(v.get_variant("chrA", 2).get_allele_string(1) == "GAAGCCAGTTCCCCGAGACGGCCAAA");
	REQUIRE(v.get_variant("chrA", 2).get_allele_string(2) == "GAAGCCAGTTCCCCTACGGCCAAA");
	REQUIRE(v.get_variant("chrA", 2).nr_of_paths() == 5);

	REQUIRE(v.get_variant("chrA", 3).get_allele_string(0) == "ACGTCCGTTCAGCCTTAGC");
	REQUIRE(v.get_variant("chrA", 3).get_allele_string(1) == "ACGTCCGTTTAGCCTTAGC");

	REQUIRE(v.get_variant("chrA", 4).get_allele_string(0) == "CCGATTTTCTTGTGCTATA");
	REQUIRE(v.get_variant("chrA", 4).get_allele_string(1) == "CCGATTTTCCTGTGCTATA");

	REQUIRE(v.get_variant("chrA", 5).get_allele_string(0) == "GGAGGGTATGAAGCCATCAC");
	REQUIRE(v.get_variant("chrA", 5).get_allele_string(1) == "GGAGGGTATTCAGCCATCAC");

	REQUIRE(v.get_variant("chrA", 6).get_allele_string(0) == "TGTGGACTTATTTGGCTAA");
	REQUIRE(v.get_variant("chrA", 6).get_allele_string(1) == "TGTGGACTTGTTTGGCTAA");

	REQUIRE(v.get_variant("chrB", 0).get_allele_string(0) == "CCACTTCATCAAGACACAA");
	REQUIRE(v.get_variant("chrB", 1).get_allele_string(0) == "GAGTATTTTGATCATAAAT");
	// there should be 4+1 paths (since the reference path is added)
	REQUIRE(v.nr_of_paths() == 5);

	v.write_path_segments("small1-segments.fa");
}

TEST_CASE("VariantReader get_overhang", "[VariantReader get_overhang]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	VariantReader v(vcf, fasta, 10, false);

	DnaSequence left_overhang;
	DnaSequence right_overhang;
	v.get_left_overhang("chrA", 0, 20, left_overhang);
	v.get_right_overhang("chrA", 0, 20, right_overhang);
	REQUIRE(left_overhang.to_string() == "TTTGGTGATCTGGAATTCCG");
	REQUIRE(right_overhang.to_string() == "CATAAGTTATGCTAAAAAAT");

	v.get_left_overhang("chrA", 1, 20, left_overhang);
	v.get_right_overhang("chrA", 1, 20, right_overhang);
	REQUIRE(left_overhang.to_string() == "GTCTGTTAAGACCTTAGCTA");
	REQUIRE(right_overhang.to_string() == "GAAGCCAGT");

	v.get_left_overhang("chrA", 2, 20, left_overhang);
	v.get_right_overhang("chrA", 2, 20, right_overhang);
	REQUIRE(left_overhang.to_string() == "GAAGCCAGT");
	REQUIRE(right_overhang.to_string() == "ACGGCCAAAACATACCATTT");
	REQUIRE(v.nr_of_paths() == 4);
}

TEST_CASE("VariantReader write_path_segments", "[VariantReader write_path_segments]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF file
	VariantReader v(vcf, fasta, 10, false);
	v.write_path_segments("../tests/data/small1-segments.fa");

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

TEST_CASE("VariantReader write_path_segments_no_variants", "[VariantReader write_path_segments]") {
	string vcf = "../tests/data/empty.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF
	VariantReader v(vcf, fasta, 10, false);
	v.write_path_segments("../tests/data/empty-segments.fa");
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

TEST_CASE("VariantReader write_genotypes_of", "[VariantReader write_genotypes_of]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";

	// read variants from VCF file
	VariantReader v(vcf, fasta, 10, false, "HG0");

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

	v.open_genotyping_outfile("../tests/data/small1-genotypes.vcf");
	v.write_genotypes_of("chrA", genotypes_chrA);
	v.write_genotypes_of("chrB", genotypes_chrB);
	v.close_genotyping_outfile();

	v.open_phasing_outfile("../tests/data/small1-phasing.vcf");
	v.write_phasing_of("chrA", genotypes_chrA);
	v.write_phasing_of("chrB", genotypes_chrB);
	v.close_genotyping_outfile();
}

TEST_CASE("VariantReader broken_vcfs", "[VariantReader broken_vcfs]") {
	string no_paths = "../tests/data/no-paths.vcf";
	string malformatted = "../tests/data/malformatted-vcf1.vcf";
	string fasta = "../tests/data/small1.fa";

	CHECK_THROWS(VariantReader(no_paths, fasta, 10, false));
	CHECK_THROWS(VariantReader(malformatted, fasta, 10, false));
}

TEST_CASE("VariantReader no-alt-alleles", "[VariantReader no-alt-alleles]") {
	string vcf = "../tests/data/no-alt-alleles.vcf";
	string fasta = "../tests/data/small1.fa";

	VariantReader v (vcf, fasta, 10, false);
	// should have skipped variant for which no alt alleles are given
	REQUIRE(v.size_of("chrA") == 1);
}

TEST_CASE("VariantReader overlapping variants", "[VariantReader overlapping variants]") {
	string vcf = "../tests/data/overlapping-variants.vcf";
	string fasta = "../tests/data/small1.fa";

	VariantReader v (vcf, fasta, 10, false);
	// should have skipped variant that is contained in another
	REQUIRE(v.size_of("chrA") == 1);
}

TEST_CASE("VariantReader get_chromosomes", "[VariantReader get_chromosomes]") {
	string vcf1 = "../tests/data/small1.vcf";
	string vcf2 = "../tests/data/small2.vcf";
	string fasta = "../tests/data/small1.fa";

	VariantReader v1(vcf1, fasta, 10, false);
	vector<string> chromosomes;
	v1.get_chromosomes(&chromosomes);
	vector<string> expected1 = {"chrA", "chrB"};
	REQUIRE(chromosomes == expected1);

	VariantReader v2(vcf2, fasta, 10, false);
	chromosomes.clear();
	v2.get_chromosomes(&chromosomes);
	vector<string> expected2 = {"chrB", "chrC", "chrA"};
	REQUIRE(chromosomes == expected2);
}

TEST_CASE("VariantReader construct_index", "[VariantReader construct_index]") {
	vector<string> sequences = {"TTTTT", "AATAGTAAAGTTATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	vector<DnaSequence> alleles;
	for (auto s : sequences) {
		alleles.push_back(DnaSequence(s));
	}
	vector<unsigned char> expected = {1,0,2,3};
	REQUIRE(construct_index(alleles, true) == expected);
}

TEST_CASE("VariantReader variant_ids1", "[VariantReader variant_ids1]") {
	string vcf = "../tests/data/small1-ids.vcf";
	string fasta = "../tests/data/small1.fa";
	VariantReader v(vcf, fasta, 10, true);

	vector<string> sequences_ref = {"TGGG", "AATAGTAAAGTTATA", "GTAGATAGATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	vector<string> sequences = {"AATAGTAAAGTTATA", "GTAGATAGATA", "AATAGTAAAGTGATA", "GGGTG", "TTG"};
	map<string,string> sequence_to_id = {{"AATAGTAAAGTTATA", "var1"}, {"GTAGATAGATA", "var2"}, {"AATAGTAAAGTGATA", "var3"}, {"GGGTG", "var4"}, {"TTG", "var5:var6"}};

	vector<DnaSequence> alleles;
	for (auto s : sequences_ref) {
		alleles.push_back(DnaSequence(s));
	}
	string chromosome = "chr1";
	vector<string> variant_ids = {"var1", "var2", "var3", "var4", "var5:var6"};
	v.insert_ids(chromosome, alleles, variant_ids, true);

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
		string result = v.get_ids(chromosome, sequences, 0 , false);
		REQUIRE(result == expected);
	}
}

TEST_CASE("VariantReader variant_ids2", "[VariantReader variant_ids2]") {
	string vcf = "../tests/data/small1-ids.vcf";
	string fasta = "../tests/data/small1.fa";
	VariantReader v(vcf, fasta, 10, true);
	vector<GenotypingResult> genotypes(2);
	v.open_genotyping_outfile("../tests/data/small1-ids-genotypes.vcf");
	
	v.write_genotypes_of("chrA", genotypes);
}

TEST_CASE("VariantReader close_to_start", "[VariantReader close_to_start]") {
	string vcf = "../tests/data/close.vcf";
	string fasta = "../tests/data/close.fa";
	vector<GenotypingResult> genotypes(1);
	VariantReader v(vcf, fasta, 31, true);
	v.open_genotyping_outfile("../tests/data/small1-ids-close.vcf");
	v.write_genotypes_of("chr10", genotypes);
}


TEST_CASE("VariantReader non_existing_path", "[VariantReader non_existing_path]") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	VariantReader v(vcf, fasta, 10, false);
	CHECK_THROWS(v.write_path_segments("nonexistent/paths_segments.fasta"));
}


TEST_CASE("VariantReader too_large_panel", "[VariantReader too_large_panel]") {
	string vcf = "../tests/data/large-panel.vcf";
	string fasta = "../tests/data/small1.fa";
	// there are more than 256 paths in the VCF, the implementation cannot handle this and should throw an error
	CHECK_THROWS(VariantReader (vcf, fasta, 10, false));
	CHECK_THROWS(VariantReader (vcf, fasta, 10, true));
}


TEST_CASE("VariantReader too_many_alleles", "[VariantReader too_many_alleles]") {
	string vcf = "../tests/data/many-alleles.vcf";
	string fasta = "../tests/data/small1.fa";
	// there are more than 256 alleles in the VCF, the implementation cannot handle this and should throw an error
	CHECK_THROWS(VariantReader (vcf, fasta, 10, false));
}


