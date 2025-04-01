#ifndef KMERPATH16_HPP
#define KMERPATH16_HPP

#include <vector>
#include <stdint.h>
#include  <string>
#include <iostream>
#include <cereal/access.hpp>
#include <cereal/types/vector.hpp>

/** Represents a sequence of kmers. **/

class KmerPath16 {
public:
	KmerPath16();
	/** indicate presence of kmer at index **/
	void set_position(unsigned short index);
	/** check given position **/
	unsigned int get_position(unsigned short index) const;
	/** compute number of kmers on this path **/
	size_t nr_kmers() const;
	friend std::ostream& operator<< (std::ostream& stream, const KmerPath16& cna);
	std::string convert_to_string() const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(offset, kmers);
	}

private:
	/** use one unsigned int to store the assignments of 32 kmers **/
	unsigned short offset;
	uint16_t kmers;
	friend cereal::access;
};

#endif // KMERPATH16_HPP
