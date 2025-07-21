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
	long double genome = 1.0L * genome_kmers;
	while (reader.next()){
		long double count = 1.0L * reader.val();
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

	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, largest_peak);

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
