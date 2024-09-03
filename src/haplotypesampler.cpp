#include "haplotypesampler.hpp"
#include <algorithm>
#include <cassert>
#include <map>
#include <math.h>

using namespace std;


HaplotypeSampler::HaplotypeSampler(vector<shared_ptr<UniqueKmers>>* unique_kmers, size_t size)
	:unique_kmers(unique_kmers)
{}

void HaplotypeSampler::get_column_minima(std::vector<unsigned int>& column, std::vector<bool>& mask, size_t& first_id, size_t& second_id, unsigned int& first_val, unsigned int& second_val)
{
 
	assert (column.size() > 1);
	assert (column.size() == mask.size());
 
    first_val = std::numeric_limits<unsigned int>::max();
	second_val = std::numeric_limits<unsigned int>::max();

    for (size_t i = 0; i < column.size(); i++) {
		// if masked with false, ignore.
		if (!mask[i]) continue;

        /* If current element is smaller than first
        then update both first and second */
        if (column[i] < first_val) {
            second_val = first_val;
			second_id = first_id;
			first_val = column[i];
            first_id = i;
        }
 
        /* If arr[i] is in between first and second
        then update second */
        else if (column[i] < second_val && i != first_id)
            second_val = column[i];
			second_id = i;
    }
}


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

void HaplotypeSampler::compute_viterbi_path() {
	size_t column_count = this->unique_kmers->size();
	init(this->viterbi_columns, column_count);
	init(this->viterbi_backtrace_columns, column_count);

	// perform Viterbi algorithm
	size_t k = (size_t) sqrt(column_count);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {
		compute_viterbi_column(column_index);
		// store sparse table. Check if previous column needs to be deleted.
		if ((k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->viterbi_columns[column_index-1];
			this->viterbi_columns[column_index-1] = nullptr;
			delete this->viterbi_backtrace_columns[column_index-1];
			this->viterbi_backtrace_columns[column_index-1] = nullptr;
		}
	}

	// find the best value in the last column
	size_t best_index = 0;
	long double best_value = 0.0L;
	DPColumn* last_column = this->viterbi_columns.at(column_count-1);
	assert (last_column != nullptr);
	for (size_t i = 0; i < last_column->column.size(); ++i) {
		long double entry = last_column->column.at(i);
		if (entry >= best_value) {
			best_value = entry;
			best_index = i;
		}
	}

	// backtracking
	this->sampled_paths.sampled_paths.push_back(vector<size_t>(column_count));
	size_t column_index = column_count - 1;
	while(true) {
		// columns might have to be re-computed
		if (this->viterbi_backtrace_columns[column_index] == nullptr) {
			size_t j = column_index / k*k;
			assert (this->viterbi_columns[j] != nullptr);
			for (j = j+1; j <= column_index; ++j) {
				compute_viterbi_column(j);
			}
		}
		// store the best path
		this->sampled_paths.sampled_paths[this->sampled_paths.sampled_paths.size()-1][column_index] = (best_index);

		if (column_index == 0) break;

		// update the best index
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);
		column_index -= 1;
	}
}

void HaplotypeSampler::compute_viterbi_column(size_t column_index) {
	// TODO implement this
	assert (column_index < this->unique_kmers->size());

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column
	DPColumn* previous_column = nullptr;

	// number of paths
	size_t nr_paths = this->unique_kmers->at(column_index)->get_nr_paths();

	// TODO: class that computes Transition probabilities?

	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index - 1];
	}

	DPColumn* current_column = new DPColumn();

	// TODO: class for computing emission probabilities?

	// backtrace column
	vector<size_t>* backtrace_column = new vector<size_t>();

	// determine indices of not yet covered haplotypes (boolean vector)
	vector<bool> current_indexes = this->sampled_paths.mask_indexes(column_index, nr_paths);

	// TODO: precompute minima for each index in current column
	// TODO: fill DP column based on precomputed minima



}

void HaplotypeSampler::update_unique_kmers() {}
