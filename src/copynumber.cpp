#include <stdexcept>
#include <sstream>
#include <functional>
#include <algorithm>
#include "copynumber.hpp"

using namespace std;

CopyNumber::CopyNumber()
:probabilities({1.0,0.0,0.0})
{}

CopyNumber::CopyNumber(long double cn_0, long double cn_1, long double cn_2, long double scaling_factor)
	:probabilities(3,0.0L)
{
	this->probabilities[0] = cn_0 / scaling_factor;
	this->probabilities[1] = cn_1 / scaling_factor;
	this->probabilities[2] = cn_2 / scaling_factor;
}

long double CopyNumber::get_probability_of(int cn){
	if( (cn < 0) || (cn > 2) ){
		ostringstream oss;
		oss << "CopyNumber::get_probability_of: Invalid copy number: " << cn;
		throw runtime_error(oss.str());
	}
	return probabilities[cn];
}

bool CopyNumber::operator==(const CopyNumber &other) const{
	for (size_t i = 0; i < 3; ++i){
		if (this->probabilities[i] != other.probabilities[i]){
			return false;
		}
	}
	return true;
}

bool CopyNumber::operator!=(const CopyNumber &other) const {
	return !(*this == other);
}
