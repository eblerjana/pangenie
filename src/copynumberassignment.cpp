#include <cassert>
#include <math.h>
#include "copynumberassignment.hpp"

using namespace std;

CopyNumberAssignment::CopyNumberAssignment()
	:kmers({0})
{}

CopyNumberAssignment::CopyNumberAssignment(vector<unsigned int> kmers)
	:kmers(kmers)
{}

void CopyNumberAssignment::set_position(size_t index, unsigned int value){
	assert (value < 3);
	size_t max_index = (this->kmers.size())*20 -1;
	while (index > max_index){
		this->kmers.push_back(0);
		max_index = (this->kmers.size())*20 -1;
	}
	// determine which block we have to modify
	size_t block = index / 20;
	unsigned int number = this->kmers[block];
	// determine which position we have to modify
	size_t position = index - (block * 20);
	// change the position to given value
	unsigned int factor = pow(3, position);
	unsigned int old_value = (number / factor) % 3;
	unsigned int tmp = number - (old_value * factor);
	this->kmers[block] = tmp + value * factor;
}

unsigned int CopyNumberAssignment::get_position(size_t index) const {
	size_t max_index = (this->kmers.size())*20 -1;
	if (index > max_index) {
		return 0;
	}
	size_t block = index / 20;
	unsigned int number = this->kmers[block];
	size_t position = index - (block * 20);
	unsigned int factor = pow(3, position);
	return (number / factor) % 3;
}

string CopyNumberAssignment::convert_to_string() const {
	string result = "";
	for (size_t i = 0; i < this->kmers.size(); ++i){
		unsigned int factor = 1;
		unsigned int assignment = this->kmers[i];
		for (size_t j = 0; j < 20; ++j){
			result += to_string((assignment / factor) % 3);
			factor *= 3;
		}
	}
	return result;
}

ostream& operator<< (ostream& stream, const CopyNumberAssignment& cna){
	stream << cna.convert_to_string();
	return stream;
}
