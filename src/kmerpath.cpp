#include <math.h>
#include "kmerpath.hpp"

using namespace std;

KmerPath::KmerPath()
	:kmers({})
{}

void KmerPath::set_position(size_t index){
	int max_index = int((this->kmers.size())*32) - 1;
	while ( (int) index > max_index){
		this->kmers.push_back(0);
		max_index = (this->kmers.size())*32 - 1;
	}
	// determine which block we have to modify
	size_t block = index / 32;
	// determine which position we have to modify
	size_t position = index - (block * 32);
	// change the position to given value
	this->kmers[block] |= (1 << position);
}

unsigned int KmerPath::get_position(size_t index) const {
	int max_index = int((this->kmers.size())*32) - 1;
	if ( (int) index > max_index) {
		return 0;
	}
	size_t block = index / 32;
	unsigned int number = this->kmers[block];
	size_t position = index - (block * 32);
	if ((number) & (1<<(position))) {
		return 1;
	} else {
		return 0;
	}
}

size_t KmerPath::nr_kmers() const {
	size_t result = 0;
	for (size_t i = 0; i < this->kmers.size(); ++i) {
		unsigned int assignment = this->kmers[i];
		assignment = assignment - ((assignment >> 1) & 0x55555555);
		assignment = (assignment & 0x33333333) + ((assignment >> 2) & 0x33333333);
		result += (((assignment + (assignment >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}
	return result;
}

string KmerPath::convert_to_string() const {
	string result = "";
	for (size_t i = 0; i < this->kmers.size()*32; ++i) {
		result += to_string(this->get_position(i));
	}
	return result;
}

ostream& operator<< (ostream& stream, const KmerPath& cna){
	stream << cna.convert_to_string();
	return stream;
}