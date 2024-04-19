#include "haplotypesampler.hpp"
#include <algorithm>
#include <map>


using namespace std;

HaplotypeSampler::HaplotypeSampler(vector<shared_ptr<UniqueKmers>>* unique_kmers)
	:unique_kmers(unique_kmers)
{}

// std::vector<std::shared_ptr<UniqueKmers>>*
void HaplotypeSampler::rank_haplotypes() const {
	for (auto unique_kmers : *this->unique_kmers) {
		vector<pair<unsigned int, float>> path_coverages;

		unsigned int nr_paths = unique_kmers->get_nr_paths();

		// determine fractions of present kmers per allele
		map<unsigned char, float> fractions = unique_kmers->covered_kmers_on_alleles();
		map<unsigned char, int> nr_kmers = unique_kmers->kmers_on_alleles();

		// determine fractions for each path
		for (unsigned int p = 0; p < nr_paths; ++p) {
			path_coverages.push_back(make_pair(p, fractions.at(unique_kmers->get_allele(p))));
		}

		// sort paths by coverage
		sort(path_coverages.begin(), path_coverages.end(), [](pair<unsigned int, float> &left, pair<unsigned int, float> &right) {return left.second < right.second;});

		// print out for now
		cout << "Ranked haplotypes for bubble position: " << unique_kmers->get_variant_position() << endl;
		for (auto p : path_coverages) {
			cout << "path" << p.first << "\t" << p.second << "\t" << (unsigned int) unique_kmers->get_allele(p.first) << "\t" << nr_kmers[unique_kmers->get_allele(p.first)] << endl;
		}
	}
}