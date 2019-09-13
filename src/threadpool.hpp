#ifndef THREADPOOL_HPP
#define THREADPOOL_HPP

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <list>

/** Code take from: http://www.mathematik.uni-ulm.de/numerik/pp/ss17/pp-folien-2017-05-18.pdf **/

class ThreadPool {
public:
	using Job = std::function<void()>;
	ThreadPool (size_t nr_threads);
	~ThreadPool ();
	void submit(Job job);
private:
	size_t nr_threads;
	bool finished;
	std::vector<std::thread> threads;
	std::mutex m;
	std::condition_variable cv;
	std::list<Job> jobs;
	void process_jobs ();
};

#endif //THREADPOOL_HPP
