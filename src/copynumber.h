#ifndef COPYNUMBER_H
#define COPYNUMBER_H

#include <vector>

class CopyNumber {
public:
	CopyNumber();
	CopyNumber(double cn_0, double cn_1, double cn_2);
	int get_copynumber();
	double get_probability_of(int cn);
private:
	std::vector<double> probabilities;
};
#endif // COPYNUMBER_H
