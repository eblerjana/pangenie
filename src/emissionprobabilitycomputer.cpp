#include "emissionprobabilitycomputer.h"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(size_t nr_paths)
	:current_index(0),
	 nr_paths(nr_paths)
{}

void EmissionProbabilityComputer::insert_kmer(string kmer, CopyNumber cn, std::vector<int> paths){
	size_t index = this->current_index;
	this->kmer_to_index.insert(pair<string,size_t>(kmer, index));
	this->kmer_to_copynumber.insert(pair<size_t,CopyNumber>(index,cn));
	this->kmer_to_path.insert(pair<size_t,vector<bool>>(index, vector<bool>(this->nr_paths, false)));
	for (size_t i = 0; i < paths.size(); ++i){
		this->kmer_to_path[index][paths[i]] = true;
	}
	this->current_index += 1;
}

long double EmissionProbabilityComputer::get_emission_probability(int path_id1, int path_id2){
	long double result = 1.0;
	for (map<string,size_t>::iterator it = this->kmer_to_index.begin(); it != this->kmer_to_index.end(); ++it){
		size_t index = it->second;
		bool in_p1 = this->kmer_to_path[index][path_id1];
		bool in_p2 = this->kmer_to_path[index][path_id2];
		CopyNumber cn = this->kmer_to_copynumber[index];
		if (in_p1 && in_p2){
			// kmer present in both paths
			result *= cn.get_probability_of(2);
		} else if (in_p1 || in_p2){
			// kmer present in one path
			result *= cn.get_probability_of(1);
		} else {
			// kmer in none of the paths
			result *= cn.get_probability_of(0);
		}
	}
	return result;
}
