#include <iostream>
#include "collapsedcolumnindexer.hpp"

using namespace std;

CollapsedColumnIndexer::CollapsedColumnIndexer (size_t variant_id)
	: variant_id(variant_id),
	 size(0)
{}

void CollapsedColumnIndexer::insert_allele (unsigned char allele) {
	this->alleles.push_back(allele);
	this->size += 1;
}

size_t CollapsedColumnIndexer::nr_alleles () const {
	return this->size;
}

unsigned char CollapsedColumnIndexer::get_allele (size_t allele_index) const {
	if (allele_index >= this->size) {
		throw runtime_error("ColumnIndexer::get_allele: index out of bounds.");
	}
	return this->alleles[allele_index];
}

pair<size_t,size_t> CollapsedColumnIndexer::get_allele_ids_at (size_t column_index) const {
	if (column_index >= this->size*this->size) {
		throw runtime_error("ColumnIndexer::get_allele_ids_at: index out of bounds.");
	}
	size_t p_id1 = column_index % this->nr_alleles();
	size_t p_id2 = (column_index / this->nr_alleles()) % this->nr_alleles();
	return pair<size_t,size_t> (p_id2, p_id1);
}
