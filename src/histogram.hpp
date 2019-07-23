#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <vector>
#include <string>

class Histogram {
public:
	Histogram(size_t max_value);
	void add_value(size_t value);
	void write_to_file(std::string filename) const;
	void smooth_histogram();
	void find_peaks(std::vector<size_t>& peak_ids, std::vector<size_t>& peak_values) const;
	friend std::ostream& operator<<(std::ostream& os, const Histogram& hist);
private:
	std::vector<size_t> histogram;
};
#endif // HISTOGRAM_HPP
