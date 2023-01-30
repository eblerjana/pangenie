#include "jellyfishreader.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <fstream>
#include "histogram.hpp"

using namespace std;

JellyfishReader::JellyfishReader (string readfile, size_t kmersize)
{
	// based on code from example: https://github.com/gmarcais/Jellyfish/blob/master/examples/query_per_sequence/query_per_sequence.cc
	this->ifs.open(readfile, ios::in|ios::binary);
	this->header = shared_ptr<jellyfish::file_header> (new jellyfish::file_header (this->ifs));
	if (kmersize != (this->header->key_len() / 2)) {
		ostringstream oss;
		oss << "JellyfishReader::JellyfishReader: given kmer size (" << kmersize  << ") does not match database kmer size (" << this->header->key_len() / 2  << ")." << endl;
		throw runtime_error(oss.str());
	}
	if (!this->header->canonical()) {
		ostringstream oss;
		oss << "JellyfishReader::JellyfishReader: in order to use the .jf file, kmers must have been counted using the -C option." << endl;
		throw runtime_error(oss.str());
	}
	if(!this->ifs.good()){
		ostringstream oss;
		oss << "JellyfishReader::JellyfishReader: Failed to parse header of file '" << readfile << "'" << endl;
		throw runtime_error(oss.str());
	}
	if(this->header->format() == binary_dumper::format) {
		this->binary_map = shared_ptr<jellyfish::mapped_file>( new jellyfish::mapped_file (readfile.c_str()));
		this->db = shared_ptr<binary_query> (new binary_query(this->binary_map->base() + this->header->offset(), this->header->key_len(), this->header->counter_len(), this->header->matrix(),
                    this->header->size() - 1, this->binary_map->length() - this->header->offset()));
		this->filename = readfile;
	} else {
		ostringstream oss;
		oss << "JellyfishReader::JellyfishReader: Unsupported format '" << this->header->format() << endl;
		throw runtime_error(oss.str());
	}
}

size_t JellyfishReader::getKmerAbundance(string kmer){
	jellyfish::mer_dna jelly_kmer(kmer);
	jelly_kmer.canonicalize();
	return this->db->check(jelly_kmer);
}

size_t JellyfishReader::getKmerAbundance(jellyfish::mer_dna jelly_kmer){
	jelly_kmer.canonicalize();
	return this->db->check(jelly_kmer);
}

size_t JellyfishReader::computeKmerCoverage(size_t genome_kmers) {
	binary_reader reader (this->ifs, this->header.get());

	long double result = 0.0L;
	while (reader.next()){
		long double count = 1.0L * reader.val();
		long double genome = 1.0L * genome_kmers;
		result += (count/genome);
	}

	// reset ifs
	this->ifs.clear();
	this->ifs.seekg(0, ios::beg);
	return 0;
}

size_t JellyfishReader::computeHistogram(size_t max_count, bool largest_peak, string filename) {
	jellyfish::mer_dna::k(this->header->key_len() / 2);
	binary_reader reader (this->ifs, this->header.get());

	Histogram histogram(max_count);
	while (reader.next()) {
		histogram.add_value(reader.val());
	}

	// write histogram values to file
	if (filename != "") {
		histogram.write_to_file(filename);
	}
	// smooth the histogram
	histogram.smooth_histogram();
	// find peaks
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histogram.find_peaks(peak_ids, peak_values);

	// identify the largest and second largest (if it exists)
	if (peak_ids.size() == 0) {
		throw runtime_error("JellyfishReader::computeHistogram: no peak found in kmer-count histogram.");
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
	// add expected abundance counts to end of hist file
	if (filename != "") {
		ofstream histofile;
		histofile.open(filename, ios::app);
		if (!histofile.good()) {
			stringstream ss;
			ss << "JellyfishReader::computeHistogram: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
			throw runtime_error(ss.str());
		}
		histofile << "parameters\t" << kmer_coverage_estimate/2.0 << '\t' << kmer_coverage_estimate << endl;
		histofile.close();
	}

	// reset ifs
	this->ifs.clear();
	this->ifs.seekg(0, ios::beg);
	return kmer_coverage_estimate;
}

JellyfishReader::~JellyfishReader() {
	if (this->ifs.is_open()) {
		this->ifs.close();
	}
}
