#include "catch.hpp"
#include "utils.hpp"
#include "../src/graphbuilder.hpp"
#include "../src/graph.hpp"
#include "../src/fastareader.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/uniquekmercomputer.hpp"
#include "../src/stepwiseuniquekmercomputer.hpp"
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cereal/archives/binary.hpp>

using namespace std;

TEST_CASE("UniqueKmerComputer test1", "[UniqueKmerComputer test1]") {
	string reference_file = "../tests/data/small1.fa";
	string read_file = "../tests/data/reads.fa";
	string variants_file = "../tests/data/small1.vcf";
	string segments_file = "../tests/data/UniqueKmerComputerTest.fa";
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder builder(variants_file, reference_file, graph, segments_file, 31, true);
	JellyfishCounter graph_counts(segments_file, 31, 1, 3000);
	shared_ptr<JellyfishCounter> read_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(read_file, 31, 1, 3000));
	UniqueKmerComputer u(&graph_counts, read_counts, graph["chrA"], 10);

	ProbabilityTable probabilities(0, 40, 50, 1.0);
	vector<shared_ptr<UniqueKmers>> result;

	u.compute_unique_kmers(&result, &probabilities);
	
	REQUIRE(result.size() == 6);

	// make sure number of total kmers equals sum of kmers on each allele,
	// i.e. a kmer does not occur on more than one allele
 	for (auto u : result) {
		size_t total_kmers = u->size();
		size_t summed_kmers = 0;
		map<unsigned short, int> counts = u->kmers_on_alleles();
		for (auto c : counts) {
			if (c.second > 0) summed_kmers += c.second;
		}
		REQUIRE(total_kmers == summed_kmers);
	}
}

TEST_CASE("UniqueKmerComputer test2", "[UniqueKmerComputer test2]") {
	string reference_file = "../tests/data/small1.fa";
	string read_file = "../tests/data/reads.fa";
	string variants_file = "../tests/data/small5.vcf";
	string segments_file = "../tests/data/UniqueKmerComputerTest.fa";
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder builder(variants_file, reference_file, graph, segments_file, 31, true);
	JellyfishCounter graph_counts(segments_file, 31, 1, 3000);
	shared_ptr<JellyfishCounter> read_counts = shared_ptr<JellyfishCounter>(new JellyfishCounter(read_file, 31, 1, 3000));
	UniqueKmerComputer u(&graph_counts, read_counts, graph["chrB"], 10);

	ProbabilityTable probabilities(0, 40, 50, 1.0);
	vector<shared_ptr<UniqueKmers>> result;

	u.compute_unique_kmers(&result, &probabilities);
	
	REQUIRE(result.size() == 1);

	// make sure number of total kmers equals sum of kmers on each allele,
	// i.e. a kmer does not occur on more than one allele
 	for (auto u : result) {
		size_t total_kmers = u->size();
		size_t summed_kmers = 0;
		map<unsigned short, int> counts = u->kmers_on_alleles();
		for (auto c : counts) {
			if (c.second > 0) summed_kmers += c.second;
		}
		REQUIRE(total_kmers == summed_kmers);
		REQUIRE(total_kmers <= 301);
	}
}


TEST_CASE("StepwiseUniqueKmerComputer test1", "[StepwiseUniqueKmerComputer test1]") {
	string reference_file = "../tests/data/small1.fa";
	string variants_file = "../tests/data/small1.vcf";
	string segments_file = "../tests/data/UniqueKmerComputerTest.fa";
	string kmers_file = "../tests/data/kmers.tsv.gz";
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder builder(variants_file, reference_file, graph, segments_file, 31, true);
	JellyfishCounter graph_counts(segments_file, 31, 1, 3000);
	StepwiseUniqueKmerComputer u(&graph_counts, graph["chrA"]);

	ProbabilityTable probabilities(0, 40, 50, 1.0);
	vector<shared_ptr<UniqueKmers>> result;

	u.compute_unique_kmers(&result, kmers_file);
	
	REQUIRE(result.size() == 6);

	// make sure number of total kmers equals sum of kmers on each allele,
	// i.e. a kmer does not occur on more than one allele
 	for (auto u : result) {
		size_t total_kmers = u->size();
		size_t summed_kmers = 0;
		map<unsigned short, int> counts = u->kmers_on_alleles();
		for (auto c : counts) {
			if (c.second > 0) summed_kmers += c.second;
		}
		REQUIRE(total_kmers == summed_kmers);
	}
}

TEST_CASE("StepwiseUniqueKmerComputer test2", "[StepwiseUniqueKmerComputer test2]") {
	string reference_file = "../tests/data/small1.fa";
	string variants_file = "../tests/data/small5.vcf";
	string segments_file = "../tests/data/UniqueKmerComputerTest.fa";
	string kmers_file = "../tests/data/kmers.tsv.gz";
	map<string, shared_ptr<Graph>> graph;
	GraphBuilder builder(variants_file, reference_file, graph, segments_file, 31, true);
	JellyfishCounter graph_counts(segments_file, 31, 1, 3000);
	StepwiseUniqueKmerComputer u(&graph_counts, graph["chrB"]);

	ProbabilityTable probabilities(0, 40, 50, 1.0);
	vector<shared_ptr<UniqueKmers>> result;

	u.compute_unique_kmers(&result, kmers_file);
	
	REQUIRE(result.size() == 1);

	// make sure number of alleles does not exceed maximum
 	for (auto u : result) {
		size_t total_kmers = u->size();
		size_t summed_kmers = 0;
		map<unsigned short, int> counts = u->kmers_on_alleles();
		for (auto c : counts) {
			if (c.second > 0) summed_kmers += c.second;
		}
		REQUIRE(total_kmers == summed_kmers);
		REQUIRE(total_kmers <= 301);
	}
}