#ifndef KMERPATH_HPP
#define KMERPATH_HPP

#include <vector>
#include "copynumberassignment.hpp"

/** Represents a sequence of kmers. **/

class KmerPath: public CopyNumberAssignment {
public:
	KmerPath();
	/** indicate presence of kmer at index **/
	void set_position(size_t index);
	/** compute number of kmers on this path **/
	size_t nr_kmers() const;
	friend CopyNumberAssignment operator+(KmerPath& p1, KmerPath& p2);
};

#endif // KMERPATH_HPP
