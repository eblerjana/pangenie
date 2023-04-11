#include <vector>
#include <string>
#include <memory>
#include "kmercounter.hpp"

void parse(std::vector<std::string>& result, std::string line, char sep);

void parse_kmer_line(std::string line, std::string& chromosome, size_t& start, std::vector<std::string>& kmers, std::vector<std::string>& flanking_kmers, bool& is_header);

unsigned short compute_local_coverage(std::vector<std::string>& kmers, std::shared_ptr<KmerCounter> read_counts, size_t kmer_coverage);