#ifndef SEQUENCE_UTILS_HPP
#define SEQUENCE_UTILS_HPP

#include <vector>
#include <cstddef>

unsigned char encode (char base);

unsigned char complement (unsigned char base);

char decode (unsigned char number);

size_t compute_kmer_coverage(std::vector<size_t>& peak_ids, std::vector<size_t>& peak_values, bool largest_peak);

#endif // SEQUENCE_UTILS_HPP
