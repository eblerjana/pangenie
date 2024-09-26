#include "haplotypesampler.hpp"
#include "samplingemissions.hpp"
#include "samplingtransitions.hpp"
#include <algorithm>
#include <cassert>
#include <map>
#include <math.h>

using namespace std;

void print_dpcolumn(DPColumn* column) {
	cout << "Print column:" << endl;
	for (size_t i = 0; i < column->column.size(); ++i) {
		cout << column->column.at(i) << endl;
	}
	cout << "--------" << endl;
}

HaplotypeSampler::HaplotypeSampler(vector<shared_ptr<UniqueKmers>>* unique_kmers, size_t size, double recombrate, long double effective_N, vector<unsigned int>* best_scores)
	:unique_kmers(unique_kmers),
	 recombrate(recombrate),
	 effective_N(effective_N)
{

	if (size < 1) return;

	// generate size Viterbi paths
	for (size_t i = 0; i < size; ++i) {
		compute_viterbi_path(best_scores);
	}

	// Update unique_kmers
	update_unique_kmers();

	// clean up
	init(this->viterbi_columns, 0);
	init(this->viterbi_backtrace_columns, 0);
	
}

void HaplotypeSampler::get_column_minima(std::vector<unsigned int>& column, std::vector<bool>& mask, size_t& first_id, size_t& second_id, unsigned int& first_val, unsigned int& second_val) const
{
 
	assert (column.size() > 1);
	assert (column.size() == mask.size());
 
    first_val = std::numeric_limits<unsigned int>::max();
	second_val = std::numeric_limits<unsigned int>::max();

	first_id = std::numeric_limits<unsigned int>::max();
	second_id = std::numeric_limits<unsigned int>::max();

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
        } else if ( (column[i] < second_val) && (i != first_id)) {
            second_val = column[i];
			second_id = i;
    	}
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

void HaplotypeSampler::compute_viterbi_path(vector<unsigned int>* best_scores) {
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
	DPColumn* last_column = this->viterbi_columns.at(column_count-1);
	size_t best_index = 0;
	unsigned int best_value = last_column->column.at(0);
	assert (last_column != nullptr);
	for (size_t i = 1; i < last_column->column.size(); ++i) {
		unsigned int entry = last_column->column.at(i);
		if (entry < best_value) {
			best_value = entry;
			best_index = i;
		}
	}

	// keep track of best DP score
	if (best_scores != nullptr) {
		best_scores->push_back(best_value);
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
	vector<bool> prev_mask;

	// number of paths
	size_t nr_paths = this->unique_kmers->at(column_index)->get_nr_paths();

	// TODO: class that computes Transition probabilities?

	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index - 1];

		// bitvector marking which path ids shall be ignored (removed in previous DP iterations)
		prev_mask = this->sampled_paths.mask_indexes(column_index-1, nr_paths-1);
	}

	DPColumn* current_column = new DPColumn();
	current_column->column = vector<unsigned int> (nr_paths);

	// TODO: class for computing emission probabilities?

	// backtrace column
	vector<size_t>* backtrace_column = new vector<size_t>(nr_paths);

	// determine indices of not yet covered haplotypes (boolean vector)
	vector<bool> current_indexes = this->sampled_paths.mask_indexes(column_index, nr_paths);

	// precompute minima for each index in current column. helper[i] contains the value of
	// the minimum value of all positions except i in previous columns.
	vector<unsigned int> helper_val(nr_paths);
	vector<unsigned int> helper_id(nr_paths);

	// currently masked indexes (removed in previous DP iterations)
	vector<bool> cur_mask = this->sampled_paths.mask_indexes(column_index, nr_paths-1);

	SamplingTransitions* transition_cost_computer = nullptr;
	SamplingEmissions emission_cost_computer(this->unique_kmers->at(column_index));


	if (column_index > 0) {
		// set up SamplingTransitions
		size_t from_variant = this->unique_kmers->at(column_index-1)->get_variant_position();
		size_t to_variant = this->unique_kmers->at(column_index)->get_variant_position();
		transition_cost_computer = new SamplingTransitions(from_variant, to_variant, this->recombrate, nr_paths, this->effective_N);

		// compute smallest and second smallest element in previous column
		size_t first_id, second_id;
		unsigned int first_val, second_val;
		this->get_column_minima(previous_column->column, prev_mask, first_id, second_id, first_val, second_val);
		// fill helper vector
		for (size_t i = 0; i < nr_paths; ++i) {
			if (cur_mask[i]) {
				if (i == first_id) {
					// need to consider second smallest value from previous column
					helper_val[i] = second_val;
					helper_id[i] = second_id;
				} else {
					// need to consider minimum from previous column
					helper_val[i] = first_val;
					helper_val[i] = first_id;
				}
			} else {
				helper_val[i] = numeric_limits<unsigned int>::max();
				helper_id[i] = numeric_limits<unsigned int>::max();
			}
		}
	}

	// TODO: check for overflows (see WH code)!!

	// fill DP column based on precomputed minima
	for (size_t i = 0; i < nr_paths; ++i) {
		// if current index is masked, skip.
		if (!cur_mask[i]) {
			current_column->column[i] = numeric_limits<unsigned int>::max();
			continue;
		}
		unsigned int previous_cell = 0;
		if (column_index > 0) {
			// check of previous value exists for same path (might be masked)
			// keep track of where the minimum came from and store in backtrace table
			previous_cell = helper_val[i] +transition_cost_computer->compute_transition_cost(true);

			// check if there was an overflow
			if (previous_cell < helper_val[i]) previous_cell = numeric_limits<unsigned int>::max();

			backtrace_column->operator[](i) = helper_id[i];

			if (prev_mask[i]) {
				unsigned int same = previous_column->column.at(i) + transition_cost_computer->compute_transition_cost(false);

				// check if there was an overflow
				if (same < previous_column->column.at(i)) same = numeric_limits<unsigned int>::max();

				if (same < previous_cell) {
					previous_cell = same;
					backtrace_column->operator[](i) = i;
				}
			}
		}
		// add Emission costs
		size_t allele = this->unique_kmers->at(column_index)->get_allele(i);
		current_column->column[i] = previous_cell + emission_cost_computer.get_emission_cost(allele);

		// check if there was an overflow
		if (current_column->column[i] < previous_cell) current_column->column[i] = numeric_limits<unsigned int>::max();
	}

	// store the column and clean up
	this->viterbi_columns.at(column_index) = current_column;
	this->viterbi_backtrace_columns.at(column_index) = backtrace_column;

	if (transition_cost_computer != nullptr) {
		delete transition_cost_computer;
	}
}

void HaplotypeSampler::update_unique_kmers() {
	// change the UniqueKmers objects so that they contain only the paths that are left after sampling
	size_t nr_paths = this->sampled_paths.sampled_paths.size();
	size_t nr_columns = this->unique_kmers->size();
	for (size_t i = 0; i < nr_columns; ++i) {
		vector<unsigned short> p(nr_paths);
		// iterate all paths sampled at this position
		for (size_t j = 0; j < nr_paths; ++j) {
			// insert sampled paths
			p[j] = this->sampled_paths.sampled_paths[j][i];
		}
		cout << "Before update" << endl;
		cout << *this->unique_kmers->operator[](i) << endl;
		this->unique_kmers->operator[](i)->update_paths(p);
		cout << "After update" << endl;
		cout << *this->unique_kmers->operator[](i) << endl;
		cout << endl;
	}
}


SampledPaths HaplotypeSampler::get_sampled_paths() const {
	return this->sampled_paths;
}