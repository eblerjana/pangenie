#include "catch.hpp"
#include "utils.hpp"
#include "../src/recombinationmap.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("RecombinationMap read_genetic_map", "[RecombinationMap read_genetic_map]"){

	string filename = "../tests/data/recomb_map.map";
	vector<pair<size_t, float>> g;
	load_genetic_map(filename, &g);

	for (auto s : g) {
		cout << s.first << " " << s.second << endl;
	}
}