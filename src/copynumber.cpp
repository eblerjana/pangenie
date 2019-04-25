#include <stdexcept>
#include <sstream>
#include "copynumber.h"

using namespace std;

CopyNumber::CopyNumber()
:probabilities({1.0,0.0,0.0})
{}

CopyNumber::CopyNumber(double cn_0, double cn_1, double cn_2)
	:probabilities({cn_0, cn_1, cn_2})
{
	// determine actual copynumber
	vector<double> tmp = {cn_0, cn_1, cn_2};
	int best_index = 0;
	double best = 0.0;
	for (size_t i = 0; i < 3; ++i){
		if (tmp[i] > best){
			best = tmp[i];
			best_index = i;
		}
	}
	probabilities.push_back(best_index);
}

int CopyNumber::get_copynumber(){
	return probabilities[3];
}

double CopyNumber::get_probability_of(int cn){
	if( (cn < 0) || (cn > 2) ){
		ostringstream oss;
		oss << "Invalid copy number: " << cn;
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
