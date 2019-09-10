#include "jellyfishcounter.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <fstream>
#include "histogram.hpp"

using namespace std;

JellyfishCounter::JellyfishCounter (string readfile, size_t kmer_size, size_t nr_threads)
{
	jellyfish::mer_dna::k(kmer_size); // Set length of mers
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = nr_threads; // TODO Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool canonical = true; // Use canonical representation

	// create the hash
	this->jellyfish_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

	// convert the string to char**
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

	// count kmers
	mer_counter jellyfish_counter(num_threads, (*jellyfish_hash), &args[0], (&args[0])+1, canonical);
	jellyfish_counter.exec_join(num_threads);

	// delete the char**
	for(size_t i = 0; i < args.size(); i++)
		delete[] args[i];

}

size_t JellyfishCounter::getKmerAbundance(string kmer){

	jellyfish::mer_dna jelly_kmer(kmer);
	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	jf_ary->get_val_for_key(jelly_kmer, &val);

	const auto end = jf_ary->end();
	return val;
}

size_t JellyfishCounter::getKmerAbundance(jellyfish::mer_dna jelly_kmer){

	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	jf_ary->get_val_for_key(jelly_kmer, &val);

	const auto end = jf_ary->end();
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

size_t JellyfishCounter::computeHistogram(size_t max_count, string filename) {
	Histogram histogram(max_count);
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto end = jf_ary->end();
	for (auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		histogram.add_value(key_val.second);
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
	// identify the two largest peaks and return rightmost one
	if (peak_ids.size() < 2) {
		throw runtime_error("JellyfishCounter: less than 2 peaks found.");
	}
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
	// add expected abundance counts to end of hist file
	if (filename != "") {
		ofstream histofile;
		histofile.open(filename, ios::app);
		histofile << "parameters\t" << 0.9 << '\t' << second_id/2.0 << '\t' << second_id << endl;
		histofile.close();
	}
	return second_id;
}

JellyfishCounter::~JellyfishCounter() {
	delete this->jellyfish_hash;
	this->jellyfish_hash = nullptr;
}
