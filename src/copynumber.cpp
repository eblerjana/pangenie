#include <stdexcept>
#include <sstream>
#include <functional>
#include <algorithm>
#include "copynumber.hpp"
#include <iostream>

using namespace std;

CopyNumber::CopyNumber()
:probabilities({1.0,0.0})
{}

CopyNumber::CopyNumber(long double cn_0, long double cn_1, long double cn_2)
	:probabilities(3,0.0L)
{
	this->probabilities[0] = cn_0;
	this->probabilities[1] = cn_1;
	this->probabilities[2] = cn_2;
}

CopyNumber::CopyNumber(long double cn_0, long double cn_1, long double cn_2, long double regularization_const)
	:probabilities(2,0.0L)
{
	long double sum = cn_0 + cn_1 + cn_2 + 3.0L * regularization_const;
	this->probabilities[0] = (cn_0 + regularization_const) / sum;
	this->probabilities[1] = (cn_1 + regularization_const) / sum;
}

long double CopyNumber::get_probability_of(int cn) const {
	if( (cn < 0) || (cn > 2) ){
		ostringstream oss;
		oss << "CopyNumber::get_probability_of: Invalid copy number: " << cn;
		throw runtime_error(oss.str());
	}
	if( (cn == 2) && (this->probabilities.size() < 3)) {
		return 1.0L - probabilities.at(0) - probabilities.at(1);
	} else {
		return probabilities.at(cn);
	}
}

bool CopyNumber::operator==(const CopyNumber &other) const{
	vector<long double> p1 = this->probabilities;
	if (p1.size() < 3) p1.push_back(1.0L - p1[0] - p1[1]);
	vector<long double> p2 = other.probabilities;
	if (p2.size() < 3) p1.push_back(1.0L - p2[0] - p2[1]);
	for (size_t i = 0; i < 3; ++i){
		if (p1 != p2){
			return false;
		}
	}
	return true;
}

bool CopyNumber::operator!=(const CopyNumber &other) const {
	return !(*this == other);
}
