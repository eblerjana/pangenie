#include "jellyfishcounter.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <fstream>
#include "histogram.hpp"

using namespace std;

vector<char*> to_args(string readfile) {
	vector<char*> args;
	istringstream iss (readfile);
	string token;
	while(iss >> token) {
		char* arg = new char [token.size() + 1];
		copy(token.begin(), token.end(), arg);
		arg[token.size()] = '\0';
		args.push_back(arg);
	}
	args.push_back(0);
	return args;
}


JellyfishCounter::JellyfishCounter (string readfile, size_t kmer_size, size_t nr_threads, uint64_t hash)
{
	jellyfish::mer_dna::k(kmer_size); // Set length of mers
	const uint64_t hash_size    = hash; // Initial size of hash, default = 3000000000.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = nr_threads; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool canonical = true; // Use canonical representation

	// create the hash
	this->jellyfish_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

	// convert the readfile to char**
	vector<char*> args = to_args(readfile);

	// count kmers
	mer_counter jellyfish_counter(num_threads, (*jellyfish_hash), &args[0], (&args[0])+1, canonical, COUNT);
	jellyfish_counter.exec_join(num_threads);

	// delete the readfile char**
	for(size_t i = 0; i < args.size(); i++)
		delete[] args[i];

}

JellyfishCounter::JellyfishCounter (string readfile, vector<string> kmerfiles, size_t kmer_size, size_t nr_threads, uint64_t hash)
{
	jellyfish::mer_dna::k(kmer_size); // Set length of mers
	const uint64_t hash_size    = hash; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = nr_threads; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool canonical = true; // Use canonical representation

	// create the hash
	this->jellyfish_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

	// convert the filenames to char**
	vector<char*> reads_args = to_args(readfile);

	// process input kmers contained in the provided FASTQ files
	for (auto kmerfile : kmerfiles) {
		vector<char*> kmer_args = to_args(kmerfile);
		mer_counter jellyfish_counter(num_threads, (*jellyfish_hash), &kmer_args[0], (&kmer_args[0])+1, canonical, PRIME);
		jellyfish_counter.exec_join(num_threads);

		// delete the kmerfile char**
		for(size_t i = 0; i < kmer_args.size(); i++)
			delete[] kmer_args[i];
	}

	// process read kmers
	mer_counter jellyfish_counter(num_threads, (*jellyfish_hash), &reads_args[0], (&reads_args[0])+1, canonical, UPDATE);
	jellyfish_counter.exec_join(num_threads);

	// delete the readfile char**
	for(size_t i = 0; i < reads_args.size(); i++)
		delete[] reads_args[i];

}

size_t JellyfishCounter::getKmerAbundance(string kmer){

	jellyfish::mer_dna jelly_kmer(kmer);
	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	jf_ary->get_val_for_key(jelly_kmer, &val);
	return val;
}

size_t JellyfishCounter::getKmerAbundance(jellyfish::mer_dna jelly_kmer){

	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	jf_ary->get_val_for_key(jelly_kmer, &val);
	return val;
}

size_t JellyfishCounter::computeKmerCoverage(size_t genome_kmers) {
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto end = jf_ary->end();
	long double result = 0.0L;
	for (auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		long double count = 1.0L * key_val.second;
		long double genome = 1.0L * genome_kmers;
		result += (count/genome);
	}
	return (size_t) ceil(result);
}

size_t JellyfishCounter::computeHistogram(size_t max_count, bool largest_peak, string filename) {
	Histogram histogram(max_count);
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto end = jf_ary->end();
	for (auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		if (key_val.second > 0) histogram.add_value(key_val.second);
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
		throw runtime_error("JellyfishCounter::computeHistogram: no peak found in kmer-count histogram.");
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
			ss << "JellyfishCounter::computeHistogram: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
			throw runtime_error(ss.str());
		}
		histofile << "parameters\t" << kmer_coverage_estimate/2.0 << '\t' << kmer_coverage_estimate << endl;
		histofile.close();
	}
	return kmer_coverage_estimate;
}

JellyfishCounter::~JellyfishCounter() {
	delete this->jellyfish_hash;
	this->jellyfish_hash = nullptr;
}
