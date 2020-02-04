#include "utils.hpp"
#include <math.h>

bool doubles_equal(double a, double b) {
	return std::abs(a - b) < 0.0000001;
}

bool compare_vectors (std::vector<double>& v1, std::vector<double>& v2) {
	if (v1.size() != v2.size()) return false;

	for (size_t i = 0; i < v1.size(); ++i) {
		if (! doubles_equal(v1[i], v2[i])) return false;
	}
	return true;
}
