#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <string>
#include "commandlineparser.hpp"
#include "variantreader.hpp"
#include "timer.hpp"
#include "threadpool.hpp"
#include "recombinationmap.hpp"

using namespace std;

struct GeneticMaps {
	mutex genmap_mutex;
	map<string, shared_ptr<RecombinationMap>> maps;
	map<string, double> runtimes;
};

int main(int argc, char* argv[]){
    CommandLineParser argument_parser;
    argument_parser.add_command("RecombOutput -v <variants.vcf> -r <reference.fasta> -m <maps.tsv>");
    argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
    argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('m', "TSV file specifying paths to chromosome-specific genetic maps (in map format). File format: <chromosome> <path/to/genetic.map>");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	string vcffile = argument_parser.get_argument('v');
	string reffile = argument_parser.get_argument('r');
    string genetic_maps = argument_parser.get_argument('m');

    // read VCF
	VariantReader variant_reader (vcffile, reffile, 31, true, "sample");
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);

	// prepare recombination maps
	GeneticMaps recombination_maps;

	map<string, string> map_filenames;
	parse_map_inputs(genetic_maps, map_filenames);

    for (auto chromosome : chromosomes) {
		if (map_filenames.find(chromosome) == map_filenames.end()) {
            cerr << "Constructing uniform recombination map for chromosome " << chromosome << endl;
			recombination_maps.maps[chromosome] = shared_ptr<RecombinationMap> (new RecombinationMap(1.26));
		} else {
            cerr << "Reading recombination map for chromosome " << chromosome << " from file: " << map_filenames[chromosome] << endl;
            recombination_maps.maps[chromosome] = shared_ptr<RecombinationMap> (new RecombinationMap(map_filenames[chromosome], &variant_reader, chromosome));
        }
    }

    // print recombination probs
    cout << "chromosome\tposition\tuniform_prob\tmap_prob" << endl;
    RecombinationMap uniform(1.26);

    for (auto chromosome : chromosomes) {
        for (size_t v  = 0; v < variant_reader.size_of(chromosome)-1; ++v) {
            size_t position = variant_reader.get_variant(chromosome, v).get_start_position();
            size_t next_position = variant_reader.get_variant(chromosome, v+1).get_start_position();
            cout << chromosome << "\t" << position << "\t" << uniform.compute_recombination_probability(position, v, next_position, v+1) << "\t" << recombination_maps.maps[chromosome]->compute_recombination_probability(position, v, next_position, v+1) << endl;
        }
    }
}