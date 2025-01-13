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
	vector<unsigned short> path_to_allele = {0, 0};
	shared_ptr<UniqueKmers> u1 =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned short> a1 = {0};
	vector<unsigned short> a2 = {1};

	path_to_allele = {1, 0};
	shared_ptr<UniqueKmers> u2 = shared_ptr<UniqueKmers>(new UniqueKmers (3000, path_to_allele));
	u2->set_undefined_allele(0);
	REQUIRE (u2->is_undefined_allele(0));
	u2->insert_kmer(20, a2);
	u2->insert_kmer(1, a2);

	// check fractions computed by UniqueKmers object
	map<unsigned short, float> fractions;
	fractions[0] = u1->fraction_present_kmers_on_allele(0);
	REQUIRE(fractions.size() == 1);
	REQUIRE(doubles_equal(fractions[0], 1.0));

	fractions.clear();
	fractions[0] = u2->fraction_present_kmers_on_allele(0);
	fractions[1] = u2->fraction_present_kmers_on_allele(1);
	REQUIRE(fractions.size() == 2);
	REQUIRE(doubles_equal(fractions[0], 1.0));
	REQUIRE(doubles_equal(fractions[1], 0.5));

	SamplingEmissions s1 (u1);
	REQUIRE(s1.get_emission_cost(0) == 0);

	SamplingEmissions s2 (u2);
	REQUIRE(s2.get_emission_cost(0) == 50);
	REQUIRE(s2.get_emission_cost(1) == 3);
}

TEST_CASE("SamplingEmissions get_emission_cost2", "[SamplingEmissions get_emission_cost2]")
{
	vector<unsigned short> path_to_allele = {0, 1};
	shared_ptr<UniqueKmers> u1 =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned short> a1 = {0};
	vector<unsigned short> a2 = {1};
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
	map<unsigned short, float> fractions;
	fractions[0] = u1->fraction_present_kmers_on_allele(0);
	fractions[1] = u1->fraction_present_kmers_on_allele(1);

	REQUIRE(doubles_equal(fractions[0], 2.0/3.0));
	REQUIRE(doubles_equal(fractions[1], 1.0));

	fractions.clear();
	fractions[0] = u2->fraction_present_kmers_on_allele(0);
	fractions[1] = u2->fraction_present_kmers_on_allele(1);

	REQUIRE(fractions.size() == 2);
	REQUIRE(doubles_equal(fractions[0], 0.0));
	REQUIRE(doubles_equal(fractions[1], 1.0/3.0));

	SamplingEmissions s1 (u1);
	REQUIRE(s1.get_emission_cost(0) == 1);
	REQUIRE(s1.get_emission_cost(1) == 0);

	SamplingEmissions s2 (u2);
	REQUIRE(s2.get_emission_cost(0) == 50);
	REQUIRE(s2.get_emission_cost(1) == 4);
}


TEST_CASE("SamplingEmissions get_emission_cost3", "[SamplingEmissions get_emission_cost3]")
{
	vector<unsigned short> path_to_allele = {0, 1};
	shared_ptr<UniqueKmers> u1 =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned short> a1 = {0};
	vector<unsigned short> a2 = {1};
	u1->insert_kmer(20, a1);
	u1->insert_kmer(1, a2);


	// check fractions computed by UniqueKmers object
	map<unsigned short, float> fractions;
	fractions[0] = u1->fraction_present_kmers_on_allele(0);
	fractions[1] = u1->fraction_present_kmers_on_allele(1);

	REQUIRE(doubles_equal(fractions[0], 1.0));
	REQUIRE(doubles_equal(fractions[1], 0.0));

	SamplingEmissions s1 (u1);
	REQUIRE(s1.get_emission_cost(0) == 0);
	REQUIRE(s1.get_emission_cost(1) == 25);
}

TEST_CASE("SamplingEmissions undefined_allele", "[SamplingEmissions undefined_allele]") {
	vector<unsigned short> path_to_allele = {0,1,2};
	shared_ptr<UniqueKmers> u =  shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned short> a1 = {0};
	vector<unsigned short> a2 = {1};
	vector<unsigned short> a3 = {2};
	u->set_undefined_allele(1);
	u->insert_kmer(20, a1);
	u->insert_kmer(2, a3);

	// check fractions computed by UniqueKmers object
	map<unsigned short, float> fractions;
	fractions[0] = u->fraction_present_kmers_on_allele(0);
	fractions[1] = u->fraction_present_kmers_on_allele(1);
	fractions[2] = u->fraction_present_kmers_on_allele(2);

	REQUIRE(fractions.size() == 3);
	REQUIRE(doubles_equal(fractions[0], 1.0));
	REQUIRE(doubles_equal(fractions[1], 1.0));
	REQUIRE(doubles_equal(fractions[2], 0.0));

	SamplingEmissions s(u);
	REQUIRE(s.get_emission_cost(0) == 0);
	REQUIRE(s.get_emission_cost(1) == 50);
	REQUIRE(s.get_emission_cost(2) == 25);
}