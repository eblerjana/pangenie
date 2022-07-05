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


void print_column(vector<long double>* column, ColumnIndexer* indexer) {
	for (size_t i = 0; i < column->size(); ++i) {
		pair<size_t,size_t> paths = indexer->get_path_ids_at(i);
		cout << setprecision(15) << column->at(i) << " paths: " << paths.first << " " <<  paths.second << endl;
	}
	cout << "" << endl;
}


HMM::HMM(vector<UniqueKmers*>* unique_kmers, ProbabilityTable* probabilities, bool run_genotyping, bool run_phasing, double recombrate, bool uniform, long double effective_N, vector<unsigned short>* only_paths, bool normalize)
	:unique_kmers(unique_kmers),
	 probabilities(probabilities),
	 genotyping_result(unique_kmers->size()),
	 recombrate(recombrate),
	 uniform(uniform),
	 effective_N(effective_N)
{
	// index all columns with at least one alternative allele
	index_columns(only_paths);

	size_t size = this->column_indexers.size();
	// initialize forward normalization sums
	this->forward_normalization_sums = vector<long double>(size, 0.0L);
	this->previous_backward_column = nullptr;

	if (run_genotyping) {
		compute_forward_prob();
		compute_backward_prob();

		if (normalize) {
			for (size_t i = 0; i < this->genotyping_result.size(); ++i) {
				genotyping_result[i].normalize();
			}
		}
	}

	if (run_phasing) {
		compute_viterbi_path();
	}
}

HMM::~HMM(){
	init(this->forward_columns,0);
	if (this->previous_backward_column != nullptr) delete this->previous_backward_column;
	init(this->viterbi_columns,0);
	init(this->viterbi_backtrace_columns,0);
	init(this->column_indexers, 0);
}

void HMM::index_columns(vector<unsigned short>* only_paths) {
	size_t column_count = this->unique_kmers->size();
	// do one forward pass to compute ColumnIndexers
	for (size_t column_index = 0; column_index < column_count; ++ column_index) {
		// get path ids of current column
		vector<unsigned short> current_paths;
		vector<unsigned char> current_alleles;
		this->unique_kmers->at(column_index)->get_path_ids(current_paths, current_alleles, only_paths);
		unsigned short nr_paths = current_paths.size();

		if (nr_paths == 0) {
			ostringstream oss;
			oss << "HMM::index_columns: column " << column_index << " is not covered by any paths.";
			throw runtime_error(oss.str());
		}

		// check whether there are any non-reference alleles in panel
		bool all_absent = true;
		for (unsigned short i = 0; i < nr_paths; ++i) {
			if ((current_alleles[i] != 0) && (!this->unique_kmers->at(column_index)->is_undefined_allele(current_alleles[i])) ) all_absent = false;
		}

		if (!all_absent) {
			// the ColumnIndexer to be filled
			ColumnIndexer* column_indexer = new ColumnIndexer(column_index);
			for (unsigned short i = 0; i < nr_paths; ++i) {
				column_indexer->insert_path(current_paths[i], current_alleles[i]);
			}
			this->column_indexers.push_back(column_indexer);
		}
	}
}

void HMM::compute_forward_prob() {
	size_t column_count = this->column_indexers.size();
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
	size_t column_count = this->column_indexers.size();
	if (column_count == 0) return;
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
	size_t column_count = this->column_indexers.size();
	if (column_count == 0) return;
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
		pair<unsigned short, unsigned short> path_ids = this->column_indexers.at(column_index)->get_path_ids_at(best_index);
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
		size_t variant_id = this->column_indexers.at(column_index)->get_variant_id();
		this->genotyping_result.at(variant_id).add_first_haplotype_allele(allele1);
		this->genotyping_result.at(variant_id).add_second_haplotype_allele(allele2);

		if (column_index == 0) break;

		// update best index 
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);
		column_index -= 1;
	}
}

void HMM::compute_forward_column(size_t column_index) {
	assert(column_index < this->column_indexers.size());
	size_t variant_id = this->column_indexers.at(column_index)->get_variant_id();

	// check whether column was computed already
	if (this->forward_columns[column_index] != nullptr) return;

	// get previous column and previous path ids (if existent)
	vector<long double>* previous_column = nullptr;
	ColumnIndexer* previous_indexer = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	
	// get ColumnIndexer
	ColumnIndexer* column_indexer = column_indexers.at(column_index);
	assert (column_indexer != nullptr);
	// nr of paths
	unsigned short nr_paths = column_indexer->nr_paths();
	
	if (column_index > 0) {
		previous_column = this->forward_columns[column_index-1];
		previous_indexer = this->column_indexers.at(column_index-1);
		size_t prev_index = this->column_indexers.at(column_index-1)->get_variant_id();
		size_t cur_index = this->column_indexers.at(column_index)->get_variant_id();
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
		
	}

	// construct new column
	vector<long double>* current_column = new vector<long double>();

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);

	// normalization
	long double normalization_sum = 0.0L;

	// state index
	size_t i = 0;
	unsigned short nr_prev_paths = 0;
	if (column_index > 0) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			unsigned short path1 = column_indexer->get_path(path_id1);
			unsigned short path2 = column_indexer->get_path(path_id2);
			long double previous_cell = 0.0L;
			if (column_index > 0) {
				// previous state index
				size_t j = 0;
				// iterate over all pairs of previous paths
				for (unsigned short prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (unsigned short prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// forward probability of previous cell
						long double prev_forward = previous_column->at(j);
						// paths corresponding to path indices
						unsigned short prev_path1 = previous_indexer->get_path(prev_path_id1);
						unsigned short prev_path2 = previous_indexer->get_path(prev_path_id2);

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
	if (normalization_sum > 0.0L) {
		this->forward_normalization_sums.at(column_index) = normalization_sum;
	} else {
		this->forward_normalization_sums.at(column_index) = 1.0L;
	}

	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}

void HMM::compute_backward_column(size_t column_index) {
	size_t column_count = this->column_indexers.size();
	assert(column_index < column_count);
	size_t variant_id = this->column_indexers.at(column_index)->get_variant_id();

	// get previous indexers and probabilitycomputers
	ColumnIndexer* previous_indexer = nullptr;
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	EmissionProbabilityComputer* emission_probability_computer = nullptr;
	vector<long double>* forward_column = this->forward_columns.at(column_index);
	
	// get ColumnIndexer
	ColumnIndexer* column_indexer = column_indexers.at(column_index);
	assert (column_indexer != nullptr);

	// nr of paths
	unsigned short nr_paths = column_indexer->nr_paths();

	if (column_index < column_count-1) {
		assert (this->previous_backward_column != nullptr);
		size_t prev_index = this->column_indexers.at(column_index)->get_variant_id();
		size_t cur_index = this->column_indexers.at(column_index+1)->get_variant_id();
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);	
		previous_indexer = this->column_indexers.at(column_index+1);
		emission_probability_computer = new EmissionProbabilityComputer(this->unique_kmers->at(this->column_indexers.at(column_index+1)->get_variant_id()), this->probabilities);

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

	// normalization
	long double normalization_sum = 0.0L;

	// normalization of forward-backward
	long double normalization_f_b = 0.0L;

	// state index
	size_t i = 0;
	unsigned short nr_prev_paths = 0;
	if (column_index < column_count - 1) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			unsigned short path1 = column_indexer->get_path(path_id1);
			unsigned short path2 = column_indexer->get_path(path_id2);
			// get alleles on current paths
			unsigned char allele1 = column_indexer->get_allele(path_id1);
			unsigned char allele2 = column_indexer->get_allele(path_id2);
			long double current_cell = 0.0L;
			if (column_index < column_count - 1) {
				// iterate over previous column (ahead of this)
				size_t j = 0;
				for (unsigned short prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (unsigned short prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// paths corresponding to path indices
						unsigned short prev_path1 = previous_indexer->get_path(prev_path_id1);
						unsigned short prev_path2 = previous_indexer->get_path(prev_path_id2);
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
			this->genotyping_result.at(variant_id).add_to_likelihood(allele1, allele2, forward_backward_prob * this->forward_normalization_sums.at(column_index));
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

	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}

void HMM::compute_viterbi_column(size_t column_index) {
	assert(column_index < this->column_indexers.size());
	size_t variant_id = this->column_indexers.at(column_index)->get_variant_id();

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column and previous path ids (if existent)
	vector<long double>* previous_column = nullptr;
	ColumnIndexer* previous_indexer = nullptr;
	
	// get ColumnIndexer
	ColumnIndexer* column_indexer = this->column_indexers.at(column_index);
	assert (column_indexer != nullptr);
	// nr of paths
	unsigned short nr_paths = column_indexer->nr_paths();
	
	TransitionProbabilityComputer* transition_probability_computer = nullptr;
	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index-1];
		previous_indexer = this->column_indexers.at(column_index-1);
		size_t prev_index = this->column_indexers.at(column_index-1)->get_variant_id();
		size_t cur_index = this->column_indexers.at(column_index)->get_variant_id();
		size_t prev_pos = this->unique_kmers->at(prev_index)->get_variant_position();
		size_t cur_pos = this->unique_kmers->at(cur_index)->get_variant_position();
		transition_probability_computer = new TransitionProbabilityComputer(prev_pos, cur_pos, this->recombrate, nr_paths, this->uniform, this->effective_N);
	}

	// construct new column
	vector<long double>* current_column = new vector<long double>();

	// emission probability computer
	EmissionProbabilityComputer emission_probability_computer(this->unique_kmers->at(variant_id), this->probabilities);

	// normalization 
	long double normalization_sum = 0.0L;

	// backtrace table
	vector<size_t>* backtrace_column = new vector<size_t>();

	// state index
	size_t i = 0;
	unsigned short nr_prev_paths = 0;
	if (column_index > 0) nr_prev_paths = previous_indexer->nr_paths();
	// iterate over all pairs of current paths
	for (unsigned short path_id1 = 0; path_id1 < nr_paths; ++path_id1) {
		for (unsigned short path_id2 = 0; path_id2 < nr_paths; ++path_id2) {
			// get paths corresponding to path indices
			unsigned short path1 = column_indexer->get_path(path_id1);
			unsigned short path2 = column_indexer->get_path(path_id2);
			long double previous_cell = 0.0L;
			if (column_index > 0) {
				// previous state index
				size_t j = 0;
				long double max_value = 0.0L;
				size_t max_index = 0;
				// iterate over all pairs of previous paths
				for (unsigned short prev_path_id1 = 0; prev_path_id1 < nr_prev_paths; ++prev_path_id1) {
					for (unsigned short prev_path_id2 = 0; prev_path_id2 < nr_prev_paths; ++prev_path_id2) {
						// paths corresponding to path indices
						unsigned short prev_path1 = previous_indexer->get_path(prev_path_id1);
						unsigned short prev_path2 = previous_indexer->get_path(prev_path_id2);
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
	
	if (transition_probability_computer != nullptr) {
		delete transition_probability_computer;
	}
}

vector<GenotypingResult> HMM::get_genotyping_result() const {
	return this->genotyping_result;
}

vector<GenotypingResult> HMM::move_genotyping_result() {
	return move(this->genotyping_result);
}
