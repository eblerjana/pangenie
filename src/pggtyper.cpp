#include <iostream>
#include "kmercounter.hpp"
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	// TODO: parse command line to get the parameters

	clock_t clock_start = clock();
	cerr << "This is PGGTyper." << endl;
	string readfile;
	string reffile;
	string vcffile;
	size_t kmersize;
	string outname;
	string sample_name;

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, kmersize, sample_name);
	string segment_file = outname + "_path_segments.fasta";
	cerr << "Write path segments to file: " << segment_file << " ..." << endl;
	variant_reader.write_path_segments(segment_file);

	// determine total genome size
	size_t genome_kmers = variant_reader.nr_of_genomic_kmers();

	// determine chromosomes present in VCF
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);
	cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

	// determine kmer copynumbers in reads
	cerr << "Count kmers in reads ..." << endl;
	KmerCounter read_kmer_counts (readfile, kmersize);
	cerr << "Compute kmer-coverage ..." << endl;
	size_t kmer_coverage = read_kmer_counts.computeKmerCoverage(genome_kmers);
	cerr << "Computed kmer-coverage: " << kmer_coverage << endl;

	// count kmers in allele + reference sequence
	cerr << "Count kmers in genome ..." << endl;
	KmerCounter genomic_kmer_counts (segment_file, kmersize);

	// prepare output files
	variant_reader.open_genotyping_outfile(outname + "_genotyping.vcf");
	variant_reader.open_phasing_outfile(outname + "_phasing.vcf");

	for (auto& chromosome : chromosomes) {
		cerr << "Processing chromosome " << chromosome << "." << endl;
		cerr << "Determine unique kmers ..." << endl;
		// determine sets of kmers unique to each variant region
		UniqueKmerComputer kmer_computer(&genomic_kmer_counts, &read_kmer_counts, &variant_reader, chromosome, kmer_coverage);
		std::vector<UniqueKmers>* unique_kmers = kmer_computer.compute_unique_kmers();

		// TESTING
//		for (size_t i = 0; i < unique_kmers->size(); ++i) {
//			cout << unique_kmers->at(i) << endl;
//		}
		// 

		// get variants on this chromosome
		cerr << "Construct HMM" << endl;
		HMM hmm(unique_kmers, variant_reader.get_variants_on_chromosome(chromosome));

		// output the genotyping results
		cerr << "Write genotyping output ..." << endl;
		variant_reader.write_genotypes_of(chromosome, hmm.get_genotyping_result());

		// output the phasing results
		cerr << "Write phasing output ..." << endl;
		variant_reader.write_phasing_of(chromosome, hmm.get_genotyping_result());
	}

	variant_reader.close_genotyping_outfile();
	variant_reader.close_phasing_outfile();

	// total time
	double cpu_time = (double)(clock() - clock_start) / CLOCKS_PER_SEC;
	cerr << "Total CPU time: " << cpu_time << " sec" << endl;

}


