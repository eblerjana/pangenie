#include <math.h>
#include "kmerpath.hpp"

using namespace std;

KmerPath::KmerPath()
{}

void KmerPath::set_position(size_t index){
	CopyNumberAssignment::set_position(index, 1);
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
