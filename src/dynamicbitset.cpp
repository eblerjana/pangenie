#include "dynamicbitset.hpp"

using namespace std;

size_t popcount64(uint64_t x ) {
	x -= (x >> 1) & 0x5555555555555555;
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
	return (x * 0x0101010101010101) >> 56;
}

DynamicBitset::DynamicBitset() 
	:bitstrings({0})
{}

bool DynamicBitset::is_set (size_t index) const {
	size_t max_index = this->bitstrings.size()*64 - 1;
	if (index < max_index) {
		size_t chunk = index / 64;
		size_t position_in_bitstring = index - (64*chunk);
		return this->bitstrings.at(chunk) & (1<<position_in_bitstring);
	} else {
		return false;
	}
}

void DynamicBitset::set (size_t index) {
	this->set_index_to(index, 1);
}


void DynamicBitset::set_index_to (size_t index, unsigned int bit) {
	size_t max_index = this->bitstrings.size()*64 - 1;
	while (index > max_index) {
		this->bitstrings.push_back(0);
		max_index = this->bitstrings.size()*64 - 1;
	}
	size_t chunk = index / 64;
	size_t position_in_bitstring = index - (64*chunk);
	if (bit == 1) {
		this->bitstrings[chunk] |= (1<<position_in_bitstring);
	} else {
		this->bitstrings[chunk] &= ~(1<<position_in_bitstring);
	}
}

void DynamicBitset::unset (size_t index) {
	this->set_index_to (index, 0);
}

void DynamicBitset::unset (size_t index, bool& was_set) {
	was_set = this->is_set(index);
	this->set_index_to (index, 0);
}

size_t DynamicBitset::count () const {
	size_t result = 0;
	for (auto b : this->bitstrings) {
		result += popcount64(b);
	}
	return result;
}

DynamicBitset operator&(DynamicBitset& a, DynamicBitset& b) {
	// perform & componentwise
	vector<uint64_t> combined;
	if (a.bitstrings.size() > b.bitstrings.size()) {
		combined = a.bitstrings;
		for (size_t i = 0; i < b.bitstrings.size(); ++i) {
			combined[i] = a.bitstrings[i] & b.bitstrings[i];
		}
	} else {
		combined = b.bitstrings;
		for (size_t i = 0; i < a.bitstrings.size(); ++i) {
			combined[i] = a.bitstrings[i] & b.bitstrings[i];
		}
	}
	DynamicBitset result;
	result.bitstrings = combined;
	return result;
}

string DynamicBitset::convert_to_string() const {
	string result = "";
	for (size_t i = 0; i < this->bitstrings.size(); ++i) {
		uint64_t paths = this->bitstrings.at(i);
		for (size_t j = 0; j < 64; ++j) {
			result += to_string(paths & 1);
			paths = paths >> 1;
		}
	}
	return result;
}

ostream& operator<< (ostream& stream, const DynamicBitset& b) {
	stream << b.convert_to_string();
	return stream;
}
