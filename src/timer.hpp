#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

class Timer {
public:
	Timer();
	double get_total_time() const;
	double get_interval_time();
private:
	std::chrono::steady_clock::time_point start_clock;
	std::chrono::steady_clock::time_point last_stop;
};

#endif //TIMER_HPP
