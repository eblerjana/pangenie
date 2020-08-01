#include <iostream>
#include "columnindexer.hpp"

using namespace std;

ColumnIndexer::ColumnIndexer (size_t variant_id)
	:variant_id(variant_id)
{}

void ColumnIndexer::insert_path (unsigned short path, unsigned char allele) {
	this->paths.push_back(path);
	this->alleles.push_back(allele);
}

unsigned short ColumnIndexer::nr_paths () const {
	return this->paths.size();
}

unsigned short ColumnIndexer::get_path (unsigned short path_index) const {
	if (path_index >= this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_allele: index out of bounds.");
	}
	return this->paths[path_index];
}

unsigned char ColumnIndexer::get_allele (unsigned short path_index) const {
	if (path_index >= this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_allele: index out of bounds.");
	}
	return this->alleles[path_index];
}

pair<unsigned short,unsigned short> ColumnIndexer::get_path_ids_at (size_t column_index) const {
	if (column_index >= this->paths.size()*this->paths.size()) {
		throw runtime_error("ColumnIndexer::get_path_ids_at: index out of bounds.");
	}
	size_t p_id1 = column_index % this->nr_paths();
	size_t p_id2 = (column_index / this->nr_paths()) % this->nr_paths();
	return pair<unsigned short,unsigned short> (p_id2, p_id1);
}

size_t ColumnIndexer::get_variant_id() const {
	return this->variant_id;
}
