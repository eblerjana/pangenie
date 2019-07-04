#include <iostream>
#include "columnindexer.hpp"

using namespace std;

ColumnIndexer::ColumnIndexer(size_t variant_id, size_t column_size)
	:variant_id(variant_id),
	 size(0),
	 index_to_paths(column_size)
{}

void ColumnIndexer::insert(pair<size_t,size_t> p) {
	if (size == this->index_to_paths.size()) {
		throw runtime_error("ColumnIndexer::insert: no more elements can be inserted.");
	}
	this->index_to_paths[this->size] = p;
	this->size += 1;
}

pair<size_t, size_t> ColumnIndexer::get_paths (size_t index) {
	if (index >= this->size) {
		throw runtime_error("ColumnIndexer::get_paths: index out of bounds.");
	}
	return this->index_to_paths[index];
}

size_t ColumnIndexer::get_size() {
	return this->size;
}
