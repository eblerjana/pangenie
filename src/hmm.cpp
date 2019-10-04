#include <utility>
#include <math.h>
#include <cassert>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "hmm.hpp"
#include "emissionprobabilitycomputer.hpp"

#include <iostream>

using namespace std;

/**
void print_column(vector<long double>* column, ColumnIndexer* indexer) {
	for (size_t i = 0; i < column->size(); ++i) {
		pair<size_t,size_t> paths = indexer->get_paths(i);
		cout << setprecision(15) << column->at(i) << " paths: " << paths.first << " " <<  paths.second << endl;
	}
	cout << "" << endl;
}
**/

HMM::HMM(vector<UniqueKmers*>* unique_kmers, bool run_genotyping, bool run_phasing, double recombrate, bool uniform, long double effective_N)
	:unique_kmers(unique_kmers),
	 genotyping_result(unique_kmers->size())
{
	size_t size = this->unique_kmers->size();

	// construct TransitionProbabilityComputers
	init(this->transition_prob_computers, size-1);
	for (size_t i = 1; i < size; ++i) {
		size_t prev_pos = this->unique_kmers->at(i-1)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(i)->get_variant_position();
		size_t nr_paths = this->unique_kmers->at(i)->get_nr_paths();
		TransitionProbabilityComputer* t = new TransitionProbabilityComputer(prev_pos, cur_pos, recombrate, nr_paths, uniform, effective_N);
		this->transition_prob_computers.at(i-1) = t;
	}
	this->previous_backward_column = nullptr;
//	cerr << "Indexing the columns ..." << endl;
	index_columns();
	if (run_genotyping) {
//		cerr << "Computing forward probabilities ..." << endl;
		compute_forward_prob();
//		cerr << "Computing backward probabilities ..." << endl;
		compute_backward_prob();
	}

	if (run_phasing) {
//		cerr << "Computing Viterbi path ..." << endl;
		compute_viterbi_path();
	}
}

HMM::~HMM(){
	init(this->forward_columns,0);
	if (this->previous_backward_column != nullptr) delete this->previous_backward_column;
	init(this->viterbi_columns,0);
	init(this->transition_prob_computers,0);
	init(this->viterbi_backtrace_columns,0);
	init(this->column_indexers, 0);
}

void HMM::index_columns() {
	size_t column_count = this->unique_kmers->size();
	init(column_indexers, column_count);
	// do one forward pass to compute ColumnIndexers
	for (size_t column_index = 0; column_index < column_count; ++ column_index) {
		// get path ids of current column
		vector<size_t> current_paths;
		vector<unsigned char> current_alleles;
		this->unique_kmers->at(column_index)->get_path_ids(current_paths, current_alleles);
		size_t nr_paths = current_paths.size();

		if (nr_paths == 0) {
			ostringstream oss;
			oss << "HMM::index_columns: column " << column_index << " is not covered by any paths.";
			throw runtime_error(oss.str());
		}

		// the ColumnIndexer to be filled
		ColumnIndexer* column_indexer = new ColumnIndexer(column_index);
		for (size_t i = 0; i < nr_paths; ++i) {
			column_indexer->insert_path(current_paths[i], current_alleles[i]);
		}

		// store the ColummIndexer
		this->column_indexers.at(column_index) = column_indexer;
	}
}

void HMM::compute_forward_prob() {
	size_t column_count = this->unique_kmers->size();
	init(this->forward_columns, column_count);
	
	// forward pass
	size_t k = (size_t) sqrt(column_count);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {;
		compute_forward_column(column_index);
		// sparse table: check whether to delete previous column
		if ( (k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->forward_columns[column_index-1];
			this->forward_columns[column_index-1] = nullptr;
		}
	}
}

void HMM::compute_backward_prob() {
	size_t column_count = this->unique_kmers->size();
	if (this->previous_backward_column != nullptr) {
		delete this->previous_backward_column;
		this->previous_backward_column = nullptr;
	}

	// backward pass
	for (int column_index = column_count-1; column_index >= 0; --column_index) {
		compute_backward_column(column_index);
	}
}

void HMM::compute_viterbi_path() {
	size_t column_count = this->unique_kmers->size();
	init(this->viterbi_columns, column_count);
	init(this->viterbi_backtrace_columns, column_count);

	// perform viterbi algorithm
	size_t k = (size_t) sqrt(column_count);
	for (size_t column_index = 0; column_index < column_count; ++column_index) {
		compute_viterbi_column(column_index);
		// sparse table: check whether to delete previous column
		if ((k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->viterbi_columns[column_index-1];
			this->viterbi_columns[column_index-1] = nullptr;
			delete this->viterbi_backtrace_columns[column_index-1];
			this->viterbi_backtrace_columns[column_index-1] = nullptr;
		}
	}

	// find best value (+ index) in last column
	size_t best_index = 0;
	long double best_value = 0.0L;
	vector<long double>* last_column = this->viterbi_columns.at(column_count-1);
	assert (last_column != nullptr);
	for (size_t i = 0; i < last_column->size(); ++i) {
		long double entry = last_column->at(i);
		if (entry >= best_value) {
			best_value = entry;
			best_index = i;
		}
	}

	// backtracking
	size_t column_index = column_count - 1;
	while (true) {
		pair<size_t, size_t> path_ids = this->column_indexers.at(column_index)->get_path_ids_at(best_index);
		unsigned char allele1 = this->column_indexers.at(column_index)->get_allele (path_ids.first);
		unsigned char allele2 = this->column_indexers.at(column_index)->get_allele (path_ids.second);

		// columns might have to be re-computed
		if (this->viterbi_backtrace_columns[column_index] == nullptr) {
			size_t j = column_index / k*k;
			assert (this->viterbi_columns[j] != nullptr);
			for (j = j+1; j<=column_index; ++j) {
				compute_viterbi_column(j);
			}
		}

		// store resulting haplotypes
		this->genotyping_result.at(column_index).add_first_haplotype_allele(allele1);
		this->genotyping_result.at(column_index).add_second_haplotype_allele(allele2);

		if (column_index == 0) break;

		// update best index 
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);
		column_index -= 1;
	}
}

void HMM::compute_forward_column(size_t column_index) {
	assert(column_index < this->unique_kmers->size());

	// check whether column was computed already
	if (this->forward_columns[column_index] != nullptr) return;

	// get previous column and previous path ids (if existent)
	vector<long double>* previous_column = nullptr;
	ColumnIndexer* previous_indexer = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	if (column_index > 0) {
		previous_column = this->forward_columns[column_index-1];
		previous_indexer = this->column_indexers.at(column_index-1);
		transition_probability_computer = this->transition_prob_computers.at(column_index-1);
	}

	// construct new column
	vector<long double>* current_column = new vector<long double>();

	// get ColumnIndexer
	ColumnIndexer* column_indexer = column_indexers.at(column_index);
	assert (column_indexer != nullptr);

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(column_index));

	// normalization
	long double normalization_sum = 0.0L;

	// state index
	size_t i = 0;
	// nr of paths
	size_t nr_paths = column_indexer->nr_paths();
	size_t nr_prev_paths = 0;
	if (column_index > 0) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (size_t path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (size_t path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			size_t path1 = column_indexer->get_path(path_id1);
			size_t path2 = column_indexer->get_path(path_id2);
			long double previous_cell = 0.0L;
			if (column_index > 0) {
				// previous state index
				size_t j = 0;
				// iterate over all pairs of previous paths
				for (size_t prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (size_t prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// forward probability of previous cell
						long double prev_forward = previous_column->at(j);
						// paths corresponding to path indices
						size_t prev_path1 = previous_indexer->get_path(prev_path_id1);
						size_t prev_path2 = previous_indexer->get_path(prev_path_id2);

						// determine transition probability
						long double transition_prob = transition_probability_computer->compute_transition_prob(prev_path1, prev_path2, path1, path2);
						previous_cell += prev_forward * transition_prob;
						j += 1;
					}
				}
			} else {
				previous_cell = 1.0L;
			}
			// determine alleles current paths (ids) correspond to
			unsigned char allele1 = column_indexer->get_allele(path_id1);
			unsigned char allele2 = column_indexer->get_allele(path_id2);
			// determine emission probability
			long double emission_prob = emission_probability_computer.get_emission_probability(allele1,allele2);

			// set entry of current column
			long double current_cell = previous_cell * emission_prob;
			current_column->push_back(current_cell);
			normalization_sum += current_cell;
			i += 1;
		}
	}

	if (normalization_sum > 0.0L) {
		// normalize the entries in current column to sum up to 1
		transform(current_column->begin(), current_column->end(), current_column->begin(), bind(divides<long double>(), placeholders::_1, normalization_sum));
	} else {
		long double uniform = 1.0L / (long double) current_column->size();
		transform(current_column->begin(), current_column->end(), current_column->begin(),  [uniform](long double c) -> long double {return uniform;});
//		cerr << "Underflow in Forward pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

	// store the column
	this->forward_columns.at(column_index) = current_column;
}

void HMM::compute_backward_column(size_t column_index) {
	size_t column_count = this->unique_kmers->size();
	assert(column_index < column_count);

	// get previous indexers and probabilitycomputers
	ColumnIndexer* previous_indexer = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	EmissionProbabilityComputer* emission_probability_computer = nullptr;
	vector<long double>* forward_column = this->forward_columns.at(column_index);

	if (column_index < column_count-1) {
		assert (this->previous_backward_column != nullptr);
		transition_probability_computer = this->transition_prob_computers.at(column_index);
		previous_indexer = this->column_indexers.at(column_index+1);
		emission_probability_computer = new EmissionProbabilityComputer(this->unique_kmers->at(column_index+1));

		// get forward probabilities (needed for computing posteriors
		if (forward_column == nullptr) {
			// compute index of last column stored
			size_t k = (size_t)sqrt(column_count);
			size_t next = min((size_t) ( (column_index / k) * k ), column_count-1);
			for (size_t j = next+1; j <= column_index; ++j) {
				compute_forward_column(j);
			}
		}

		forward_column = this->forward_columns.at(column_index);
		assert (forward_column != nullptr);
	}

	// construct new column
	vector<long double>* current_column = new vector<long double>();

	// get ColumnIndexer
	ColumnIndexer* column_indexer = column_indexers.at(column_index);
	assert (column_indexer != nullptr);

	// normalization
	long double normalization_sum = 0.0L;

	// normalization of forward-backward
	long double normalization_f_b = 0.0L;

	// state index
	size_t i = 0;
	// nr of paths
	size_t nr_paths = column_indexer->nr_paths();
	size_t nr_prev_paths = 0;
	if (column_index < column_count - 1) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (size_t path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (size_t path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			size_t path1 = column_indexer->get_path(path_id1);
			size_t path2 = column_indexer->get_path(path_id2);
			// get alleles on current paths
			unsigned char allele1 = column_indexer->get_allele(path_id1);
			unsigned char allele2 = column_indexer->get_allele(path_id2);
			long double current_cell = 0.0L;
			if (column_index < column_count - 1) {
				// iterate over previous column (ahead of this)
				size_t j = 0;
				for (size_t prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (size_t prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// paths corresponding to path indices
						size_t prev_path1 = previous_indexer->get_path(prev_path_id1);
						size_t prev_path2 = previous_indexer->get_path(prev_path_id2);
						// alleles on previous paths
						unsigned char prev_allele1 = previous_indexer->get_allele(prev_path_id1);
						unsigned char prev_allele2 = previous_indexer->get_allele(prev_path_id2);
						long double prev_backward = this->previous_backward_column->at(j);
						// determine transition probability
						long double transition_prob = transition_probability_computer->compute_transition_prob(path1, path2, prev_path1, prev_path2);
						current_cell += prev_backward * transition_prob * emission_probability_computer->get_emission_probability(prev_allele1, prev_allele2);
						j += 1;
					}
				}
			} else {
				current_cell = 1.0L;
			}
			// store computed backward prob in column
			current_column->push_back(current_cell);
			normalization_sum += current_cell;

			// compute forward_prob * backward_prob
			long double forward_backward_prob = forward_column->at(i) * current_cell;
			normalization_f_b += forward_backward_prob;

			// update genotype likelihood
			this->genotyping_result.at(column_index).add_to_likelihood(allele1, allele2, forward_backward_prob);
			i += 1;
		}
	}

	if (normalization_sum > 0.0L) {
		transform(current_column->begin(), current_column->end(), current_column->begin(), bind(divides<long double>(), placeholders::_1, normalization_sum));
	} else {
		long double uniform = 1.0L / (long double) current_column->size();
		transform(current_column->begin(), current_column->end(), current_column->begin(), [uniform](long double c) -> long double {return uniform;});
//		cerr << "Underflow in Backward pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

//	cout << "FORWARD COLUMN: " << endl;
//	print_column(forward_column, column_indexer);

//	cout << "BACKWARD COLUMN: "  << endl;
//	print_column(current_column, column_indexer);

	// store computed column (needed for next step)
	if (this->previous_backward_column != nullptr) {
		delete this->previous_backward_column;
		this->previous_backward_column = nullptr;
	}
	this->previous_backward_column = current_column;
	if (emission_probability_computer != nullptr) delete emission_probability_computer;

	// delete forward column as it's not needed any more
	if (this->forward_columns.at(column_index) != nullptr) {
		delete this->forward_columns.at(column_index);
		this->forward_columns.at(column_index) = nullptr;
	}

	if (normalization_f_b > 0.0L) {
		// normalize the GenotypingResults likelihoods 
		this->genotyping_result.at(column_index).divide_likelihoods_by(normalization_f_b);
	}
}

void HMM::compute_viterbi_column(size_t column_index) {
	assert(column_index < this->unique_kmers->size());

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column and previous path ids (if existent)
	vector<long double>* previous_column = nullptr;
	ColumnIndexer* previous_indexer = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index-1];
		previous_indexer = this->column_indexers.at(column_index-1);
		transition_probability_computer = this->transition_prob_computers.at(column_index-1);
	}

	// construct new column
	vector<long double>* current_column = new vector<long double>();

	// get ColumnIndexer
	ColumnIndexer* column_indexer = this->column_indexers.at(column_index);
	assert (column_indexer != nullptr);

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(column_index));

	// normalization 
	long double normalization_sum = 0.0L;

	// backtrace table
	vector<size_t>* backtrace_column = new vector<size_t>();

	// state index
	size_t i = 0;
	// nr of paths
	size_t nr_paths = column_indexer->nr_paths();
	size_t nr_prev_paths = 0;
	if (column_index > 0) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (size_t path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (size_t path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			size_t path1 = column_indexer->get_path(path_id1);
			size_t path2 = column_indexer->get_path(path_id2);
			long double previous_cell = 0.0L;
			if (column_index > 0) {
				// previous state index
				size_t j = 0;
				long double max_value = 0.0L;
				size_t max_index = 0;
				// iterate over all pairs of previous paths
				for (size_t prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (size_t prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// paths corresponding to path indices
						size_t prev_path1 = previous_indexer->get_path(prev_path_id1);
						size_t prev_path2 = previous_indexer->get_path(prev_path_id2);
						// probability of previous cell
						long double prev_prob = previous_column->at(j);
						// determine transition probability
						long double transition_prob = transition_probability_computer->compute_transition_prob(prev_path1, prev_path2, path1, path2);
						prev_prob *= transition_prob;
						if (prev_prob >= max_value) {
							max_value = prev_prob;
							max_index = j;
						}
						j += 1;
					}
				}
				previous_cell = max_value;
				backtrace_column->push_back(max_index);
			} else {
				previous_cell = 1.0L;
			}

			// determine alleles current paths (ids) correspond to
			unsigned char allele1 = column_indexer->get_allele(path_id1);
			unsigned char allele2 = column_indexer->get_allele(path_id2);
			// determine emission probability
			long double emission_prob = emission_probability_computer.get_emission_probability(allele1,allele2);
			// set entry of current column
			long double current_cell = previous_cell * emission_prob;
			current_column->push_back(current_cell);
			normalization_sum += current_cell;
			i += 1;
		}
	}

	if (normalization_sum > 0.0L) {
		// normalize the entries in current column to sum up to 1 
		transform(current_column->begin(), current_column->end(), current_column->begin(), bind(divides<long double>(), placeholders::_1, normalization_sum));
	} else {
		long double uniform = 1.0L / (long double) current_column->size();
		transform(current_column->begin(), current_column->end(), current_column->begin(),  [uniform](long double c) -> long double {return uniform;});
//		cerr << "Underflow in Viterbi pass at position: " << this->unique_kmers->at(column_index)->get_variant_position() << ". Column set to uniform." << endl;
	}

	// store the column
	this->viterbi_columns.at(column_index) = current_column;
	if (column_index > 0) assert(backtrace_column->size() == column_indexer->nr_paths()*column_indexer->nr_paths());
	this->viterbi_backtrace_columns.at(column_index) = backtrace_column;
}

vector<GenotypingResult> HMM::get_genotyping_result() const {
	return this->genotyping_result;
}
