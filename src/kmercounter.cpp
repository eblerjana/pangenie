#include "kmercounter.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>

using namespace std;


KmerCounter::KmerCounter (string readfile, std::string outname, size_t kmer_size)
{
	jellyfish::mer_dna::k(kmer_size); // Set length of mers
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = 1; // TODO Number of concurrent threads
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

size_t KmerCounter::getKmerAbundance(string kmer){

	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	jf_ary->get_val_for_key(kmer, &val);
	return val;
}

size_t KmerCounter::computeKmerCoverage(size_t genome_kmers) const {
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

KmerCounter::~KmerCounter() {
	delete this->jellyfish_hash;
	this->jellyfish_hash = nullptr;
}
