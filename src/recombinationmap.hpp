#ifndef RECOMBINATIONMAP_H
#define RECOMBINATIONMAP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream> 
#include "variantreader.hpp"

void load_genetic_map(std::string filename, std::vector<std::pair<size_t,float>>* result);

class RecombinationMap {
public:
	/**
	* @param filename name of the genetic map file
	* @param variants VariantReader with information on variant positions
	* @param chromosome which chromosome is to be considered
	**/
	RecombinationMap(std::string& filename, VariantReader* variants, std::string chromosome);
	/**
	* @param constant recombination rate (uniform recombination probabilities)
	**/
	RecombinationMap(float recomb_rate);
	/** compute recombination probability between two variant positions **/
	long double compute_recombination_probability(size_t left_variant, size_t left_variant_id, size_t right_variant, size_t right_variant_id);
private:
	bool uniform;
	long double recomb_rate;
	std::vector<long double> cumulative_distances;
	void compute_recombination_cost_map(std::vector<std::pair<size_t,float>>* genetic_map, VariantReader* variants, std::string chromosome);
};
#endif // RECOMBINATIONMAP_H 
