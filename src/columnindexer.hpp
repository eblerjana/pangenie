#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>

/** 
* Keep track of paths and alleles in a column
**/

class ColumnIndexer {
public:
	ColumnIndexer(size_t variant_id);
	/** insert a path and the allele it covers **/
	void insert_path(size_t path, unsigned char allele);
	/** number of paths inserted **/
	size_t nr_paths () const;
	/** get path at index path_id **/
	size_t get_path (size_t path_index) const;
	/** get allele at index path_id **/
	unsigned char get_allele (size_t path_index) const;
	/** get column index a pair of states corresponds to **/
	std::pair<size_t,size_t> get_path_ids_at (size_t column_index) const;

private:
	size_t variant_id;
	size_t size;
	std::vector<size_t> paths;
	std::vector<unsigned char> alleles;
};

#endif // COLUMNINDEXER_HPP
