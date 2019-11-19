#ifndef DYNAMICBITSET_H
#define DYNAMICBITSET_H

#include <vector>
#include <cstdint>
#include <iostream>
#include <string>

class DynamicBitset {
public:
	DynamicBitset();
	/**
	* @param bits bitstring represented as string object
	**/
	DynamicBitset(std::string bits);
	/** check if bit at position index is set to true **/
	bool is_set (size_t index) const;
	/** set bit at position index to true **/
	void set (size_t index);
	/** set bit at position index to false **/
	void unset (size_t index);
	/** set bit at position index to false and report if the bit previously was set to true **/
	void unset (size_t index, bool& was_set);
	/** count the number of bits set to true **/
	size_t count () const;
	/** convert bitset to string **/
	std::string convert_to_string() const;
	/** logical and operator **/
	friend DynamicBitset operator&(DynamicBitset& a, DynamicBitset& b);
	/** << operator **/
	friend std::ostream& operator<<(std::ostream& stream, const DynamicBitset& b);

private:
	std::vector<uint64_t> bitstrings;
	void set_index_to (size_t index, unsigned int bit); 
};

# endif // DYNAMICBITSET_H
