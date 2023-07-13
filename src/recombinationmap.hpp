#ifndef RECOMBINATIONMAP_H
#define RECOMBINATIONMAP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream> 
#include "variantreader.hpp"

void parse_map_inputs(std::string filename, std::map<std::string, std::string>& output_maps);
void load_genetic_map(std::string filename, std::vector<std::pair<size_t,long double>>* result);

class RecombinationMap {
public:
	/**
	* @param filename name of the genetic map file
	* @param variants VariantReader with information on variant positions
	* @param chromosome which chromosome is to be considered
	**/
	RecombinationMap(std::string& filename, VariantReader* variants, std::string chromosome);
	/**
	* @param recomb_rate constant recombination rate (uniform recombination probabilities)
	**/
	RecombinationMap(long double recomb_rate);
	/** compute recombination probability between two variant positions 
	* @param left_variant position of left variant. Only considered in case of uniform recombination rate
	* @param left_variant_id index of left variant.
	* @param right_variant position of right variant. Only considered in case on uniform recombination rate
	* @param right_variant_id index of right variant
	**/
	long double compute_recombination_probability(size_t left_variant, size_t left_variant_id, size_t right_variant, size_t right_variant_id) const;
private:
	bool uniform;
	long double recomb_rate;
	std::vector<long double> cumulative_distances;
	void compute_recombination_cost_map(std::vector<std::pair<size_t,long double>>* genetic_map, VariantReader* variants, std::string chromosome);
};
#endif // RECOMBINATIONMAP_H 
