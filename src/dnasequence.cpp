#include <stdexcept>
#include "dnasequence.hpp"

using namespace std;

unsigned char encode (char base) {
	switch (base) {
		case 'A': return 0;
		case 'a': return 0;
		case 'C': return 1;
		case 'c': return 1;
		case 'G': return 2;
		case 'g': return 2;
		case 'T': return 3;
		case 't': return 3;
		default: return 4;
	}
}

char decode (unsigned char number) {
	switch (number) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N';
	}
}

DnaSequence::DnaSequence() 
	:size(0)
{}

DnaSequence::DnaSequence(string& sequence)
	:size(0)
{
	this->append(sequence);	
}

void DnaSequence::append(string& sequence) {
	for (auto base : sequence) {
		unsigned char number = encode(base);
		if (this->size % 2 == 0) {
			this->sequence.push_back(number << 4);
		} else {
			size_t index = this->size / 2;
			unsigned char updated = (this->sequence[index]) | number;
			this->sequence[index] = updated;
		}
		this->size += 1;
	}
}

char DnaSequence::operator[](size_t position) const {
	if (position >= this->size) {
		throw runtime_error("DnaSequence::operator[]: index out of bounds.");
	}

	size_t index = position / 2;
	unsigned char number = this->sequence[index];
	if (position % 2 == 0) {
		return decode(number >> 4);
	} else {
		return decode(number & 15);
	}
}

size_t DnaSequence::length() const {
	return this->size;
}

void DnaSequence::substr(size_t start, size_t end, string& result) const {;
	result.clear();
	for (size_t i = start; i < end; ++i) {
		result += (*this)[i];
	}
}
