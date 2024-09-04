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

TEST_CASE("HaplotypeSampler get_column_minima", "[HaplotypeSampler get_column_minima]") {
	vector<unsigned int> column = {10, 2, 14, 1};
	vector<bool> mask = {true, true, true, true};

	unsigned int first_val;
	unsigned int second_val;
	size_t first_id;
	size_t second_id;

	HaplotypeSampler h (nullptr, 0);

	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 3);
	REQUIRE(second_id == 1);
	REQUIRE(first_val == 1);
	REQUIRE(second_val == 2);
}

TEST_CASE("HaplotypeSampler get_column_minima2", "[HaplotypeSampler get_column_minima2]") {
	vector<unsigned int> column = {10, 2, 14, 2};
	vector<bool> mask = {true, true, true, true};

	unsigned int first_val;
	unsigned int second_val;
	size_t first_id;
	size_t second_id;

	HaplotypeSampler h (nullptr, 0);

	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 1);
	REQUIRE(second_id == 3);
	REQUIRE(first_val == 2);
	REQUIRE(second_val == 2);
}

TEST_CASE("HaplotypeSampler get_column_minima3", "[HaplotypeSampler get_column_minima3]") {
	vector<unsigned int> column = {10, 10, 10, 10};
	vector<bool> mask = {true, true, true, true};

	unsigned int first_val;
	unsigned int second_val;
	size_t first_id;
	size_t second_id;

	HaplotypeSampler h (nullptr, 0);

	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 0);
	REQUIRE(second_id == 1);
	REQUIRE(first_val == 10);
	REQUIRE(second_val == 10);
}


TEST_CASE("HaplotypeSampler get_column_minima4", "[HaplotypeSampler get_column_minima4]") {
	vector<unsigned int> column = {10, 10, 20};
	vector<bool> mask = {true, true, true};

	unsigned int first_val;
	unsigned int second_val;
	size_t first_id;
	size_t second_id;

	HaplotypeSampler h (nullptr, 0);

	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 0);
	REQUIRE(second_id == 1);
	REQUIRE(first_val == 10);
	REQUIRE(second_val == 10);

	column = {10, 20, 20};
	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 0);
	REQUIRE(second_id == 1);
	REQUIRE(first_val == 10);
	REQUIRE(second_val == 20);
	
}

TEST_CASE("HaplotypeSampler get_column_minima5", "[HaplotypeSampler get_column_minima5]") {
	vector<unsigned int> column = {10, 20, 30};
	vector<bool> mask = {true, true, true};

	unsigned int first_val;
	unsigned int second_val;
	size_t first_id;
	size_t second_id;

	HaplotypeSampler h (nullptr, 0);
	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);


	REQUIRE(first_id == 0);
	REQUIRE(second_id == 1);
	REQUIRE(first_val == 10);
	REQUIRE(second_val == 20);

	mask = {true, false, true};
	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 0);
	REQUIRE(second_id == 2);
	REQUIRE(first_val == 10);
	REQUIRE(second_val == 30);

	mask = {false, true, true};
	h.get_column_minima(column, mask, first_id, second_id, first_val, second_val);

	REQUIRE(first_id == 1);
	REQUIRE(second_id == 2);
	REQUIRE(first_val == 20);
	REQUIRE(second_val == 30);
	
}