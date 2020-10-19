#include <chrono>
#include <future>
#include <thread>
#include <pthread.h>
#include <pybind11/pybind11.h>

extern "C" void symbolic_regress1();
extern "C" void symbolic_regress2();
extern "C" void symbolic_regress3();

template <typename TF, typename TDuration, class... TArgs>
std::result_of_t<TF&&(TArgs&&...)> runWithTimeout(TF&& f, TDuration timeout, TArgs&&... args)
{
	using R = std::result_of_t<TF&&(TArgs&&...)>;
	std::packaged_task<R(TArgs...)> task(f);
	auto future = task.get_future();
	std::thread thr(std::move(task), std::forward<TArgs>(args)...);
	auto handle = thr.native_handle();
	if (future.wait_for(timeout) != std::future_status::timeout)
	{
		thr.join();
		printf("Function exited before timeout\n");
		return future.get(); // this will propagate exception from f() if any
	}
	else
	{
		thr.detach();
		pthread_cancel(handle);
		printf("Function has been cancelled due to timeout\n");
	}
}

void sr1(unsigned int tryTimeSeconds)
{
	runWithTimeout(symbolic_regress1, std::chrono::seconds(tryTimeSeconds));
}

void sr2(unsigned int tryTimeSeconds)
{
	runWithTimeout(symbolic_regress2, std::chrono::seconds(tryTimeSeconds));
}

void sr3(unsigned int tryTimeSeconds)
{
	runWithTimeout(symbolic_regress3, std::chrono::seconds(tryTimeSeconds));
}

PYBIND11_MODULE(_feynman, m) {
    m.doc() = "Feynman symbolic regression native module";
    m.def("symbolic_regress1", &sr1);
    m.def("symbolic_regress2", &sr2);
    m.def("symbolic_regress3", &sr3);
}

