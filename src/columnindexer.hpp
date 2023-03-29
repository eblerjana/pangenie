#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>
#include <memory>
#include "uniquekmers.hpp"

/** 
* Keep track of paths and alleles in a column
**/

class ColumnIndexer {
public:
	ColumnIndexer(std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers, std::vector<unsigned short>* only_paths);
	/** get variant ID corresponding to column index **/
	size_t get_variant_id(size_t column_index) const;
	/** get the number of columns **/
	size_t size() const;
	/** number of paths **/
	unsigned short nr_paths() const;
	/** get path_id at a given index **/
	unsigned short get_path (unsigned short path_index) const;
	/** get the allele covered by a path in a column **/
	unsigned char get_allele (unsigned short path_index, size_t column_index) const;
	/** get path ids corresponding to an index inside of a column **/
	std::pair<unsigned short,unsigned short> get_path_ids_at (size_t position) const;

private:
	std::vector<size_t> variant_positions;
	std::vector<unsigned short> paths;
	std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers;
};

#endif // COLUMNINDEXER_HPP
