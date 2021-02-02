#ifndef KMERPATH_HPP
#define KMERPATH_HPP

#include <vector>
#include <stdint.h>
#include  <string>
#include <iostream>

/** Represents a sequence of kmers. **/

class KmerPath {
public:
	KmerPath();
	/** indicate presence of kmer at index **/
	void set_position(size_t index);
	/** check given position **/
	unsigned int get_position(size_t index) const;
	/** compute number of kmers on this path **/
	size_t nr_kmers() const;
	friend std::ostream& operator<< (std::ostream& stream, const KmerPath& cna);
	std::string convert_to_string() const;

private:
	/** use one unsigned int to store the assignments of 32 kmers **/
	std::vector<uint32_t> kmers;
};

#endif // KMERPATH_HPP
