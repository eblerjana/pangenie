#include "catch.hpp"
#include "utils.hpp"
#include "../src/samplingemissions.hpp"
#include "../src/uniquekmers.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace std;

TEST_CASE("SamplingEmissions get_emission_cost1", "[SamplingEmissions get_emission_cost1]")
{
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

	// check fractions computed by UniqueKmers object
	map<unsigned char, float> fractions = u1->covered_kmers_on_alleles();
	REQUIRE(fractions.size() == 1);
	REQUIRE(doubles_equal(fractions[0], 1.0));

	fractions = u2->covered_kmers_on_alleles();
	REQUIRE(fractions.size() == 2);
	REQUIRE(doubles_equal(fractions[0], 1.0));
	REQUIRE(doubles_equal(fractions[1], 0.5));

	SamplingEmissions s1 (u1);
	REQUIRE(s1.get_emission_cost(0) == 0);

	SamplingEmissions s2 (u2);
	REQUIRE(s2.get_emission_cost(0) == 0);
	REQUIRE(s2.get_emission_cost(1) == 3);
}

TEST_CASE("SamplingEmissions get_emission_cost2", "[SamplingEmissions get_emission_cost2]")
{
	vector<unsigned char> path_to_allele = {0, 1};
	shared_ptr<UniqueKmers> u1 =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1->insert_kmer(20, a1);
	u1->insert_kmer(10, a1);
	u1->insert_kmer(1, a1);
	u1->insert_kmer(3, a2);

	path_to_allele = {0, 1};
	shared_ptr<UniqueKmers> u2 = shared_ptr<UniqueKmers>(new UniqueKmers (3000, path_to_allele));
	u2->set_undefined_allele(0);
	REQUIRE (u2->is_undefined_allele(0));
	u2->insert_kmer(1, a1);
	u2->insert_kmer(1, a1);
	u2->insert_kmer(20, a2);
	u2->insert_kmer(2, a2);
	u2->insert_kmer(0, a2);

	// check fractions computed by UniqueKmers object
	map<unsigned char, float> fractions = u1->covered_kmers_on_alleles();
	REQUIRE(fractions.size() == 2);
	REQUIRE(doubles_equal(fractions[0], 2.0/3.0));
	REQUIRE(doubles_equal(fractions[1], 1.0));

	fractions = u2->covered_kmers_on_alleles();
	REQUIRE(fractions.size() == 2);
	REQUIRE(doubles_equal(fractions[0], 0.0));
	REQUIRE(doubles_equal(fractions[1], 1.0/3.0));

	SamplingEmissions s1 (u1);
	REQUIRE(s1.get_emission_cost(0) == 1);
	REQUIRE(s1.get_emission_cost(1) == 0);

	SamplingEmissions s2 (u2);
	REQUIRE(s2.get_emission_cost(0) == 25);
	REQUIRE(s2.get_emission_cost(1) == 4);
}