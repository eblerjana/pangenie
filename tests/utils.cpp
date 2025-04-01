#include "utils.hpp"
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

bool doubles_equal(double a, double b) {
	return abs(a - b) < 0.0000001;
}

bool compare_vectors (vector<double>& v1, vector<double>& v2) {
	if (v1.size() != v2.size()) return false;

	for (size_t i = 0; i < v1.size(); ++i) {
		if (! doubles_equal(v1[i], v2[i])) return false;
	}
	return true;
}

void parse_vcf_lines(string filename, vector<vector<string>>& computed_lines) {
	ifstream file(filename);
	string line;
	while(getline(file, line)) {
		vector<string> tokens;
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;
		istringstream iss(line);
		string token;
		while(getline(iss, token, '\t'))
			tokens.push_back(token);
		computed_lines.push_back(tokens);
	}
}