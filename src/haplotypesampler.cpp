#include "haplotypesampler.hpp"
#include "samplingtransitions.hpp"
#include "timer.hpp"
#include <algorithm>
#include <cassert>
#include <map>
#include <math.h>
#include <fstream>

using namespace std;

void print_dpcolumn(DPColumn* column) {
	cout << "Print column:" << endl;
	for (size_t i = 0; i < column->column.size(); ++i) {
		cout << column->column.at(i) << endl;
	}
	cout << "--------" << endl;
}

HaplotypeSampler::HaplotypeSampler(vector<shared_ptr<UniqueKmers>>* unique_kmers, size_t size, double recombrate, long double effective_N, vector<unsigned int>* best_scores, bool add_reference, string path_output, string chromosome, unsigned short allele_penalty)
	:unique_kmers(unique_kmers),
	 recombrate(recombrate),
	 effective_N(effective_N),
	 allele_penalty(allele_penalty)
{
	Timer timer;

	if (size < 1) return;

	// initialize SamplingEmissions
	this->emission_costs.reserve(unique_kmers->size());
	for (size_t column_index = 0; column_index < unique_kmers->size(); ++column_index) {
		this->emission_costs.push_back(SamplingEmissions(this->unique_kmers->at(column_index)));
	}

	cerr << "HaplotypeSampler initialize SamplingEmissions " << unique_kmers->at(0)->get_variant_position() << ": " << timer.get_interval_time() << " sec" << endl;

	// generate size Viterbi paths
	for (size_t i = 0; i < size; ++i) {
		compute_viterbi_path(best_scores);
		cerr << "HaplotypeSampler compute_viterbi_path iteration " << i << ": " << unique_kmers->at(0)->get_variant_position() << ": " << timer.get_interval_time() << " sec" << endl;
	}

	cerr << "HaplotypeSampler compute_viterbi_path total: " << unique_kmers->at(0)->get_variant_position() << ": " << timer.get_interval_time() << " sec" << endl;

	if (add_reference) this->sampled_paths.sampled_paths.push_back(vector<size_t>(unique_kmers->size(), 0));

	// print the sampled paths if requested
	if (path_output != "") {
		ofstream path_outfile;
		path_outfile.open(path_output);

		// print the header line
		path_outfile << "#chromosome\tposition";
		for (size_t path_id = 0; path_id < this->sampled_paths.sampled_paths.size(); ++path_id) {
			path_outfile << "\tHaplotypeID_path" << path_id << "\tRecombination_path" << path_id;
		}
		path_outfile << endl;

		// print stats for each position
		for (size_t column_index = 0; column_index < unique_kmers->size(); ++column_index) {
			path_outfile << chromosome << "\t" << this->unique_kmers->at(column_index)->get_variant_position();
			for (size_t path_id = 0; path_id < this->sampled_paths.sampled_paths.size(); ++path_id) {
				path_outfile << "\t" << this->sampled_paths.sampled_paths.at(path_id).at(column_index);
				path_outfile << "\t" << this->sampled_paths.recombination(column_index, path_id);
			}
			path_outfile << endl;
		}
	}

	/** 
	// DEBUGGING: write allele costs to a file
	if (debugfile != "") {
		ofstream debug_outfile;
		debug_outfile.open(debugfile);
		for (size_t column_index = 0; column_index < unique_kmers->size(); ++column_index) {
			// get transition cost
			unsigned int trans_cost = 0;
			if (column_index > 0) {
				size_t from_variant = this->unique_kmers->at(column_index-1)->get_variant_position();
				size_t to_variant = this->unique_kmers->at(column_index)->get_variant_position();
				size_t nr_paths = this->unique_kmers->at(column_index)->get_nr_paths();
				trans_cost = SamplingTransitions(from_variant, to_variant, this->recombrate, nr_paths, this->effective_N).compute_transition_cost(true);
			}

			std::map<unsigned short, int> kmers_on_alleles = this->unique_kmers->at(column_index)->kmers_on_alleles();
			bool first = true;
			debug_outfile << "variant_" << this->unique_kmers->at(column_index)->get_variant_position() << ", trans_cost: " << trans_cost << " ";
			for (auto a : kmers_on_alleles) {
				unsigned short cost = this->emission_costs.at(column_index).get_emission_cost(a.first);
				if (!first) debug_outfile << ",";
				debug_outfile << "cost:" <<  (unsigned int) cost << "(allele: " << (unsigned int) a.first  << ", #kmers:" << a.second << ")";
				first = false;
			}
			debug_outfile << endl;
		}
	}
	// DEBUGGING: END
	**/


	// Update unique_kmers
	update_unique_kmers();

	cerr << "HaplotypeSampler update_unique_kmers " << unique_kmers->at(0)->get_variant_position() << ": " << timer.get_total_time() << " sec" << endl;

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
	assert (last_column != nullptr);
	size_t best_index = 0;
	unsigned int best_value = last_column->column.at(0);
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
	vector<size_t> path(column_count);
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
		path[column_index] = best_index;

		// penalize allele covered by selected path
		unsigned short best_allele = this->unique_kmers->at(column_index)->get_allele(best_index);
		emission_costs.at(column_index).penalize(best_allele, allele_penalty);

		if (column_index == 0) break;

		// update the best index
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);

		// current column is no longer needed. Delete it.
		delete this->viterbi_columns[column_index];
		this->viterbi_columns[column_index] = nullptr;
		delete this->viterbi_backtrace_columns[column_index];
		this->viterbi_backtrace_columns[column_index] = nullptr;
		column_index -= 1;
	}
	this->sampled_paths.sampled_paths.push_back(path);
}

void HaplotypeSampler::compute_viterbi_column(size_t column_index) {
	assert (column_index < this->unique_kmers->size());

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column
	DPColumn* previous_column = nullptr;
	vector<bool> prev_mask;

	// number of paths
	size_t nr_paths = this->unique_kmers->at(column_index)->get_nr_paths();

	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index - 1];

		// bitvector marking which path ids shall be ignored (removed in previous DP iterations)
		prev_mask = this->sampled_paths.mask_indexes(column_index-1, nr_paths-1);
	}

	DPColumn* current_column = new DPColumn();
	current_column->column = vector<unsigned int> (nr_paths);

	// backtrace column
	vector<size_t>* backtrace_column = new vector<size_t>(nr_paths, numeric_limits<unsigned int>::max());

	// precompute minima for each index in current column. helper[i] contains the value of
	// the minimum value of all positions except i in previous columns.
	vector<unsigned int> helper_val(nr_paths);
	vector<unsigned int> helper_id(nr_paths);

	// currently masked indexes (removed in previous DP iterations)
	vector<bool> cur_mask = this->sampled_paths.mask_indexes(column_index, nr_paths-1);
 
	SamplingTransitions* transition_cost_computer = nullptr;

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
					helper_id[i] = first_id;
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
			previous_cell = helper_val[i] + transition_cost_computer->compute_transition_cost(true);

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
		current_column->column[i] = previous_cell + emission_costs.at(column_index).get_emission_cost(allele);

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
		this->unique_kmers->operator[](i)->update_paths(p);
	}
}


SampledPaths HaplotypeSampler::get_sampled_paths() const {
	return this->sampled_paths;
}