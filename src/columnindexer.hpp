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
	void insert_path(unsigned short path, unsigned char allele);
	/** number of paths inserted **/
	unsigned short nr_paths () const;
	/** get path at index path_id **/
	unsigned short get_path (unsigned short path_index) const;
	/** get allele at index path_id **/
	unsigned char get_allele (unsigned short path_index) const;
	/** get column index a pair of states corresponds to **/
	std::pair<unsigned short,unsigned short> get_path_ids_at (size_t column_index) const;
	/** **/
	size_t get_variant_id() const;

private:
	std::vector<unsigned short> paths;
	std::vector<unsigned char> alleles;
	size_t variant_id;
};

#endif // COLUMNINDEXER_HPP
