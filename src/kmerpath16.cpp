#include <math.h>
#include <stdexcept>
#include <cassert>
#include "kmerpath16.hpp"

using namespace std;

KmerPath16::KmerPath16()
	:offset(0),
	kmers(0)
{}

void KmerPath16::set_position(unsigned short index){
	if (this->kmers == 0) {
		// no kmers inserted yet. Determine offset.
		this->offset = index;
	}

	// check if index is valid, i.e. lies within allowed range represented by this KmerPath16 object.
	unsigned short upper_limit = this->offset + 16;
	unsigned short lower_limit = this->offset;

	if ((index < lower_limit) || (index >= upper_limit)) {
		assert(this->kmers > 0);
		throw runtime_error("KmerPath16::KmerPath16: index is invalid");
	}

	// set position
	unsigned short position = index - this->offset;
	this->kmers |= (1 << position);
}

unsigned int KmerPath16::get_position(unsigned short index) const {
	unsigned short upper_limit = this->offset + 16;
	unsigned short lower_limit = this->offset;

	if ((index < lower_limit) || (index >= upper_limit)) {
		// outside the limits, so assume absence of kmer.
		return 0;
	} else {
		unsigned short position = index - this->offset;
		if (this->kmers & (1 << position)) {
			return 1;
		} else {
			return 0;
		}
	}
}

size_t KmerPath16::nr_kmers() const {
	uint32_t assignment = this->kmers;
	assignment = assignment - ((assignment >> 1) & 0x55555555);
	assignment = (assignment & 0x33333333) + ((assignment >> 2) & 0x33333333);
	return (((assignment + (assignment >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

string KmerPath16::convert_to_string() const {
	string result = "";
	unsigned short upper_limit = this->offset + 16;
	for (size_t i = 0; i < upper_limit; ++i) {
		result += to_string(this->get_position(i));
	}
	return result;
}

ostream& operator<< (ostream& stream, const KmerPath16& cna){
	stream << cna.convert_to_string();
	return stream;
}