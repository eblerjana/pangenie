#include <vector>
#include <string>

bool doubles_equal(double a, double b);

bool compare_vectors (std::vector<double>& v1, std::vector<double>& v2);

void parse_vcf_lines(std::string filename, std::vector<std::vector<std::string>>& lines);