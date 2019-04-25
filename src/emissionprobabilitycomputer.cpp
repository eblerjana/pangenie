#include <stdexcept>
#include <sstream>
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
	this->kmer_to_path.insert(pair<size_t,vector<bool> >(index, vector<bool>(this->nr_paths, false)));
	for (size_t i = 0; i < paths.size(); ++i){
		if ( (paths[i] < 0) || (paths[i] >= this->nr_paths) ){
			ostringstream oss;
			oss << "Invalid path-ID given: " << paths[i];
			throw runtime_error(oss.str());
		}
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


bool EmissionProbabilityComputer::contains_kmer(string kmer){
	map<string,size_t>::const_iterator it = this->kmer_to_index.find(kmer);
	if (it != this->kmer_to_index.end()){
		return true;
	} else {
		return false;
	}
}

CopyNumber EmissionProbabilityComputer::get_copynumber_of(string kmer){
	if (this->contains_kmer(kmer)){
		size_t index = this->kmer_to_index[kmer];
		return this->kmer_to_copynumber[index];
	} else {
		ostringstream oss;
		oss << "Requested kmer: " << kmer << " does not exist.";
		throw runtime_error(oss.str());
	}
}

vector<int> EmissionProbabilityComputer::get_paths_of(string kmer){
	if (this->contains_kmer(kmer)){
		size_t index = this->kmer_to_index[kmer];
		vector<int> result;
		for (size_t i = 0; i < this->nr_paths; ++i){
			if (this->kmer_to_path[index][i]){
				result.push_back(i);
			}
		}
		return result;
	} else {
		ostringstream oss;
		oss << "Requested kmer: " << kmer << " does not exist.";
		throw runtime_error(oss.str());
	}
}
