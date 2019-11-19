#ifndef COLLAPSEDCOLUMNINDEXER_HPP
#define COLLAPSEDCOLUMNINDEXER_HPP

#include <utility>
#include <vector>

/** 
* Keep track of alleles in a column of the collapsed HMM
**/

class CollapsedColumnIndexer {
public:
	CollapsedColumnIndexer(size_t variant_id);
	/** insert an allele **/
	void insert_allele(unsigned char allele);
	/** number of alleles inserted **/
	size_t nr_alleles () const;
	/** get allele at index path_id **/
	unsigned char get_allele (size_t path_index) const;
	/** get column index a pair of states corresponds to **/
	std::pair<size_t,size_t> get_allele_ids_at (size_t column_index) const;

private:
	size_t variant_id;
	size_t size;
	std::vector<unsigned char> alleles;
};

#endif // COLLAPSEDCOLUMNINDEXER_HPP
