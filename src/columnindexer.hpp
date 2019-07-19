#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>

/** 
* Assign column index to pair of paths it represents.
**/

typedef std::pair< std::pair<size_t,size_t>, std::pair<unsigned char, unsigned char> > State;

class ColumnIndexer {
public:
	ColumnIndexer(size_t variant_id, size_t column_size);
	/** insert a pair of paths and corresponding alleles (state in the HMM) **/
	void insert(std::pair<size_t,size_t> p, std::pair<unsigned char, unsigned char> a);
	/** get pair of path (state) an index corresponds to **/
	std::pair<size_t, size_t> get_paths (size_t index);
	/** get pair of alleles an index corresponds to **/
	std::pair<unsigned char, unsigned char> get_alleles (size_t index);
	/** get size of the column (number of states) **/
	size_t get_size();
private:
	size_t variant_id;
	size_t size;
	/** maps index to state **/
	std::vector<State> index_to_state;
};

#endif // COLUMNINDEXER_HPP
