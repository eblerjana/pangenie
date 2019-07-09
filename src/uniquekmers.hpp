#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include "copynumber.hpp"
#include "kmerpath.hpp"

/*
* Represents the set of unique kmers for a variant position.
*/


class UniqueKmers {
public:
	UniqueKmers(size_t variant_id, size_t variant_position);
	/** insert a path that does not contain any unique kmers **/
	void insert_empty_path(size_t path);
	/** insert a kmer 
	* @param kmer sequence of kmer
	* @param cn copy number probabilities of kmer
	* @param paths Vector of paths this kmer occurs on
	**/
	void insert_kmer(CopyNumber cn, std::vector<size_t>& paths);
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	size_t get_variant_index();
	size_t get_variant_position();
	CopyNumber get_copynumber_of(size_t kmer_index);
	/** number of unique kmers **/
	size_t size() const;
	/** get all paths covering this position **/
	void get_path_ids(std::vector<size_t>& result);
	friend std::ostream& operator<< (std::ostream& stream, const UniqueKmers& uk);

private:
	size_t variant_id;
	size_t variant_pos;
	size_t current_index;
	std::vector<CopyNumber> kmer_to_copynumber;
	std::map<size_t, KmerPath> paths;
	CopyNumberAssignment combine_paths(size_t path_id1, size_t path_id2);
	friend class EmissionProbabilityComputer;
};
# endif // UNIQUEKMERS_HPP
