#include "catch.hpp"
#include "../src/haplotypesampler.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/probabilitytable.hpp"
#include "../src/probabilitycomputer.hpp"
#include "utils.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace std;

TEST_CASE("HaplotypeSampler", "[HaplotypeSampler]"){

	vector<unsigned char> path_to_allele = {0, 0};
	shared_ptr<UniqueKmers> u1 =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};

	path_to_allele = {1, 0};
	shared_ptr<UniqueKmers> u2 = shared_ptr<UniqueKmers>(new UniqueKmers (3000, path_to_allele));
	u2->set_undefined_allele(0);
	REQUIRE (u2->is_undefined_allele(0));
	u2->insert_kmer(20, a2);
	u2->insert_kmer(1, a2);

	ProbabilityTable probs (0,1,21,0.0L);
	probs.modify_probability(0, 20, CopyNumber(0.01,0.01,0.9));
	probs.modify_probability(0, 1, CopyNumber(0.9,0.3,0.1));

	vector<shared_ptr<UniqueKmers>> unique_kmers = {u1,u2};

	HaplotypeSampler haplotypesampler = HaplotypeSampler(&unique_kmers, 0);
	haplotypesampler.rank_haplotypes();

}
