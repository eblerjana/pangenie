#include <stdexcept>
#include <sstream>
#include "copynumber.hpp"

using namespace std;

CopyNumber::CopyNumber()
:probabilities({1.0,0.0,0.0})
{}

CopyNumber::CopyNumber(double cn_0, double cn_1, double cn_2)
	:probabilities({cn_0, cn_1, cn_2})
{}

double CopyNumber::get_probability_of(int cn){
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
