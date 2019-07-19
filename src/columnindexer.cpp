#include <iostream>
#include "columnindexer.hpp"

using namespace std;

ColumnIndexer::ColumnIndexer(size_t variant_id, size_t column_size)
	:variant_id(variant_id),
	 size(0),
	 index_to_state(column_size)
{}

void ColumnIndexer::insert(pair<size_t,size_t> p, pair<unsigned char, unsigned char> a) {
	if (size == this->index_to_state.size()) {
		throw runtime_error("ColumnIndexer::insert: no more elements can be inserted.");
	}
	this->index_to_state[this->size] = make_pair(p, a);
	this->size += 1;
}

pair<size_t,size_t> ColumnIndexer::get_paths (size_t index) {
	if (index >= this->size) {
		throw runtime_error("ColumnIndexer::get_paths: index out of bounds.");
	}
	return this->index_to_state.at(index).first;
}

pair<unsigned char, unsigned char> ColumnIndexer::get_alleles (size_t index) {
	if (index >= this->size) {
		throw runtime_error("ColumnIndexer::get_alleles: index out of bounds.");
	}
	return this->index_to_state.at(index).second;
}

size_t ColumnIndexer::get_size() {
	return this->size;
}
