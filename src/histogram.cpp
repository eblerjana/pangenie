#include "histogram.hpp"
#include <fstream>
#include <iostream>

using namespace std;

Histogram::Histogram(size_t max_value) 
	:histogram(max_value+1, 0)
{}

void Histogram::add_value(size_t value) {
	if (value < this->histogram.size()) {
		this->histogram[value] += 1;
	}
}

void Histogram::write_to_file(string filename) const {
	ofstream histfile;
	histfile.open(filename);
	for (size_t i = 0; i < this->histogram.size(); ++i) {
		histfile << i << '\t' << this->histogram.at(i) << endl;
	}
	histfile.close();
}

void Histogram::smooth_histogram() {
	for (size_t i = 1; i < this->histogram.size()-1; ++i) {
		this->histogram[i] = (this->histogram[i-1] + this->histogram[i] + this->histogram[i+1]) / 3;
	}
}

void Histogram::find_peaks(vector<size_t>& peak_ids, vector<size_t>& peak_values) const {
	bool direction = 0;
	size_t prev_val = 0;
	for (size_t i = 0; i < this->histogram.size(); ++i) {
		size_t value = this->histogram[i];
		if (prev_val < value) {
			direction = 0;
		} else if (prev_val > value) {
			if (direction != 1) {
				peak_ids.push_back(i-1);
				peak_values.push_back(prev_val);
			}
			direction = 1;
		}
		prev_val = value;
	}
}

ostream& operator<<(ostream& os, const Histogram& hist) {
	for (size_t i = 0; i < hist.histogram.size(); ++i) {
		os << i << '\t' << hist.histogram[i] << endl;
	}
	return os;
}
