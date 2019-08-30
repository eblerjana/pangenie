#include "timer.hpp"

using namespace std;

Timer::Timer()
	:start_clock(chrono::steady_clock::now()),
	 last_stop(chrono::steady_clock::now())
{}

double Timer::get_total_time() const {
	chrono::steady_clock::time_point now = chrono::steady_clock::now();
	return chrono::duration_cast<std::chrono::nanoseconds>(now - this->start_clock).count() / 1000000000.0;
}

double Timer::get_interval_time() {
	chrono::steady_clock::time_point now = chrono::steady_clock::now();
	double interval_time = chrono::duration_cast<std::chrono::nanoseconds>(now - this->last_stop).count();
	this->last_stop = now;
	return interval_time / 1000000000.0;
}
