#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>

/** 
* Assign column index to pair of paths it represents.
**/

class ColumnIndexer {
public:
	ColumnIndexer(size_t variant_id, size_t column_size);
	/** insert a pair of paths (state in the HMM) **/
	void insert(std::pair<size_t,size_t> p);
	/** get pair of path (state) an index corresponds to **/
	std::pair<size_t, size_t> get_paths (size_t index);
	/** get size of the column (number of states) **/
	size_t get_size();
private:
	size_t variant_id;
	size_t size;
	/** maps index to state **/
	std::vector<std::pair<size_t,size_t>> index_to_paths;
};

#endif // COLUMNINDEXER_HPP
