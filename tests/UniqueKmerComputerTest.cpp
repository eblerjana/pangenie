#include "catch.hpp"
#include "utils.hpp"
#include "../src/graphbuilder.hpp"
#include "../src/graph.hpp"
#include "../src/fastareader.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/uniquekmercomputer.hpp"
#include <string>
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

	for (auto u : result) {
		cout << *u << endl;
	}

}