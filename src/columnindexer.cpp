#include <iostream>
#include "columnindexer.hpp"

using namespace std;

ColumnIndexer::ColumnIndexer (size_t variant_id)
	: variant_id(variant_id),
	 size(0)
{}

void ColumnIndexer::insert_path (size_t path, unsigned char allele) {
	this->paths.push_back(path);
	this->alleles.push_back(allele);
	this->size += 1;
}

size_t ColumnIndexer::nr_paths () const {
	return this->size;
}

size_t ColumnIndexer::get_path (size_t path_index) const {
	if (path_index >= this->size) {
		throw runtime_error("ColumnIndexer::get_allele: index out of bounds.");
	}
	return this->paths[path_index];
}

unsigned char ColumnIndexer::get_allele (size_t path_index) const {
	if (path_index >= this->size) {
		throw runtime_error("ColumnIndexer::get_allele: index out of bounds.");
	}
	return this->alleles[path_index];
}

pair<size_t,size_t> ColumnIndexer::get_path_ids_at (size_t column_index) const {
	if (column_index >= this->size*this->size) {
		throw runtime_error("ColumnIndexer::get_path_ids_at: index out of bounds.");
	}
	size_t p_id1 = column_index % this->nr_paths();
	size_t p_id2 = (column_index / this->nr_paths()) % this->nr_paths();
	return pair<size_t,size_t> (p_id2, p_id1);
}
