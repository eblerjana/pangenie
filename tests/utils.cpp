#include "utils.hpp"
#include <math.h>

bool doubles_equal(double a, double b) {
	return std::abs(a - b) < 0.0000001;
}
