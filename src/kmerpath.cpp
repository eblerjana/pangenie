#include <math.h>
#include "kmerpath.hpp"

using namespace std;

KmerPath::KmerPath()
{}

void KmerPath::set_position(size_t index){
	CopyNumberAssignment::set_position(index, 1);
}

size_t KmerPath::nr_kmers() const {
	size_t result = 0;
	for (size_t i = 0; i < this->kmers.size(); ++i) {
		unsigned int factor = 1;
		unsigned int assignment = this->kmers[i];
		for (size_t j = 0; j < 20; ++j) {
			result += ((assignment / factor) % 3);
			factor *= 3;
		}
	}
	return result;
}

CopyNumberAssignment operator+(KmerPath& p1, KmerPath& p2){
	CopyNumberAssignment result;
	vector<unsigned int> combined;
	if (p1.kmers.size() > p2.kmers.size()){
		combined = p1.kmers;
		for (size_t i = 0; i < p2.kmers.size(); ++i){
			combined[i] += p2.kmers[i];
		}
	} else {
		combined = p2.kmers;
		for (size_t i = 0; i < p1.kmers.size(); ++i){
			combined[i] += p1.kmers[i];
		}
	}
	return CopyNumberAssignment(combined);
}
