#include <stdexcept>
#include "dnasequence.hpp"
#include <iostream>
#include "sequenceutils.hpp"

using namespace std;

DnaSequence::DnaSequence() 
	:even_length(true),
	 is_undefined(false)
{}

DnaSequence::DnaSequence(string& sequence)
	:even_length(true),
	 is_undefined(false)
{
	this->append(sequence);	
}

void DnaSequence::append(string& seq) {
	for (auto base : seq) {
		unsigned char number = encode(base);
		// check whether base was defined
		if (number == 4) this->is_undefined = true;
		if (this->even_length) {
			this->sequence.push_back(number << 4);
		} else {
			size_t index = this->size() / 2;
			unsigned char updated = (this->sequence[index]) | number;
			this->sequence[index] = updated;
		}
		this->even_length = !(this->even_length);
	}
}

void DnaSequence::append(DnaSequence seq) {
	if (this->even_length) {
		for (size_t i = 0; i < seq.sequence.size(); ++i) {
			this->sequence.push_back(seq.sequence.at(i));
		}
	} else {
		if (seq.size() == 0) return;
		unsigned char current = this->sequence.at(this->sequence.size()-1);
		this->sequence.pop_back();
		for (size_t i = 0; i < seq.size(); ++i) {
			unsigned char s = seq.sequence.at(i/2);
			if (i % 2 == 0) {
				current |= (s >> 4);
				this->sequence.push_back(current);
			} else {
				current = (s << 4);
			}
		}
		if (seq.size() % 2 == 0) this->sequence.push_back(current);
	}
	this->even_length = this->even_length == (seq.size() % 2 == 0);
	this->is_undefined = this->is_undefined || seq.contains_undefined();
}

void DnaSequence::reverse() {
	vector<unsigned char> reversed;
	if (this->even_length) {
		for (vector<unsigned char>::reverse_iterator it = this->sequence.rbegin(); it != this->sequence.rend(); ++it) {
			unsigned char reversed_element = ((*it) >> 4) | ((*it) << 4);
			reversed.push_back(reversed_element);
		}
	} else {
		unsigned char current = this->sequence[this->sequence.size() - 1];
		this->sequence.pop_back();
		for (vector<unsigned char>::reverse_iterator it = this->sequence.rbegin(); it != this->sequence.rend(); ++it) {
			unsigned char second = (*it) << 4;
			reversed.push_back( (second >> 4) | current );
			current = (*it) & 240;
		}
		reversed.push_back(current);
	}
	this->sequence = move(reversed);
}

void DnaSequence::reverse_complement() {
	vector<unsigned char> reverse_compl;
	if (this->even_length) {
		for (vector<unsigned char>::reverse_iterator it = this->sequence.rbegin(); it != this->sequence.rend(); ++it) {
			unsigned char first = complement((*it) >> 4);
			unsigned char second = complement((*it) & 15);
			reverse_compl.push_back( (second << 4) | first);
		}
	} else {
		unsigned char current = complement(this->sequence[this->sequence.size() - 1] >> 4) << 4;
		this->sequence.pop_back();
		for (vector<unsigned char>::reverse_iterator it = this->sequence.rbegin(); it != this->sequence.rend(); ++it) {
			unsigned char second = complement((*it) & 15);
			unsigned char first = complement((*it) >> 4);
			reverse_compl.push_back(current | second);
			current = (first << 4);
		}
		reverse_compl.push_back(current);
	}
	this->sequence = move(reverse_compl);
}

char DnaSequence::operator[](size_t position) const {
	if (position >= this->size()) {
		throw runtime_error("DnaSequence::operator[]: index out of bounds.");
	}

	unsigned char number = this->sequence[position / 2];
	if (position % 2 == 0) {
		return decode(number >> 4);
	} else {
		return decode(number & 15);
	}
}

DnaSequence DnaSequence::base_at(size_t position) const {
	if (position >= this->size()) {
		throw runtime_error("DnaSequence::base_at: index out of bounds.");
	}

	DnaSequence result;
	unsigned char number = this->sequence[position / 2];
	if (position % 2 == 0) {
		result.sequence.push_back(number & 240);
	} else {
		result.sequence.push_back(number << 4);
	}
	result.even_length = false;
	return result;
}

size_t DnaSequence::size() const {
	size_t result = this->sequence.size() * 2;
	if (!this->even_length) {
		result -= 1;
	}
	return result;
}

void DnaSequence::substr(size_t start, size_t end, string& result) const {
	result.clear();
	for (size_t i = start; i < end; ++i) {
		result += (*this)[i];
	}
}

void DnaSequence::substr(size_t start, size_t end, DnaSequence& result) const {
	vector<unsigned char> substring;
	if (start % 2 == 0) {
		size_t start_index = start / 2;
		size_t steps = ((end + 1) -start)/ 2;
		substring = vector<unsigned char>(this->sequence.begin() + start_index, this->sequence.begin() + start_index + steps);
		if (end % 2 != 0) substring[substring.size()-1] &= 240;
	} else {
		unsigned char current = (this->sequence[start/2]) << 4;
		for (size_t i = start+1; i < end; i+=2) {
			unsigned char element = this->sequence[i/2];
			current |= (element >> 4);
			substring.push_back(current);
			current = element <<  4;
		}
		if (end % 2 == 0) {
			substring.push_back(current);	
		}
	}

	bool undefined = false;
	if (this->is_undefined) {
		// check whether there is an undefined base in subsequence
		for (auto elem : substring) {
			if ((elem & 68) != 0) {
				undefined = true;
				break;
			}
		}
	}

	result.sequence = move(substring); 
	result.even_length = ((end - start) % 2 == 0);
	result.is_undefined = undefined;
}

string DnaSequence::to_string() const {
	string result;
	for (size_t i = 0; i < this->size(); ++i) {
		result += (*this)[i];
	}
	return result;
}

void DnaSequence::clear() {
	this->sequence.clear();
	this->even_length = true;
}

bool DnaSequence::operator<(const DnaSequence& dna) const {
	return (this->sequence < dna.sequence);
}

bool operator==(const DnaSequence& dna1, const DnaSequence& dna2) {
	return (dna1.sequence == dna2.sequence) && (dna1.size() == dna2.size());
}

bool operator!=(const DnaSequence& dna1, const DnaSequence& dna2) {
	return !(dna1 == dna2);
}

bool DnaSequence::contains_undefined() const {
	return this->is_undefined;
}
