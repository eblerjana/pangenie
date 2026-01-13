#include "sequenceutils.hpp"
#include <stdexcept>
#include <iostream>

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

unsigned char complement (unsigned char base) {
	switch(base) {
		case 0: return 3;
		case 1: return 2;
		case 2: return 1;
		case 3: return 0;
		case 4: return 4;
		default: throw runtime_error("complement: invalid base.");
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

size_t compute_kmer_coverage(vector<size_t>& peak_ids, vector<size_t>& peak_values, bool largest_peak) {
	// identify the largest and second largest (if it exists)
	if (peak_ids.size() == 0) {
		throw runtime_error("sequenceutils::computeHistogram: no peak found in kmer-count histogram.");
	}
	size_t kmer_coverage_estimate = -1;
	if (peak_ids.size() < 2) {
		cerr << "Histogram peak: " << peak_ids[0] << " (" << peak_values[0] << ")" << endl;
		kmer_coverage_estimate = peak_ids[0];
	} else {
		size_t largest, second, largest_id, second_id;
		if (peak_values[0] < peak_values[1]){
			largest = peak_values[1];
			largest_id = peak_ids[1];
			second = peak_values[0];
			second_id = peak_ids[0];
		} else {
			largest = peak_values[0];
			largest_id = peak_ids[0];
			second = peak_values[1];
			second_id = peak_ids[1];
		}
		for (size_t i = 0; i < peak_values.size(); ++i) {
			if (peak_values[i] > largest) {
				second = largest;
				second_id = largest_id;
				largest = peak_values[i];
				largest_id = peak_ids[i];
			} else if ((peak_values[i] > second) && (peak_values[i] != largest)) {
				second = peak_values[i];
				second_id = peak_ids[i];
			}
		}
		cerr << "Histogram peaks: " << largest_id << " (" << largest << "), " << second_id << " (" << second << ")" << endl;
		if (largest_peak) {
			kmer_coverage_estimate = largest_id;
		}else {
			kmer_coverage_estimate = second_id;
		}
	}
	
	return kmer_coverage_estimate;
}
