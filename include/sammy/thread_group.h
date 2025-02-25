#ifndef SAMMY_THREAD_GROUP_H_INCLUDED_
#define SAMMY_THREAD_GROUP_H_INCLUDED_

#include <cassert>
#include <condition_variable>
#include <future>
#include <mutex>
#include <thread>
#include <list>

namespace sammy {

/**
 * A group/pool of threads that handle tasks
 * which all return results of a specific type R.
 */
template <class R> class ThreadGroup {
  public:
	using ThreadIndex = std::size_t;
    using FunctionType = R(ThreadIndex);
    using TaskType = std::packaged_task<FunctionType>;
    using ResultType = R;

    explicit ThreadGroup(std::size_t desired_extra_threads) noexcept
        : m_thread_count(desired_extra_threads) {}

    explicit ThreadGroup() noexcept
        : ThreadGroup(std::thread::hardware_concurrency() - 1) {}

    // no copy or move allowed
    ThreadGroup(const ThreadGroup&) = delete;
    ThreadGroup(ThreadGroup&&) = delete;
    ThreadGroup& operator=(const ThreadGroup&) = delete;
    ThreadGroup& operator=(ThreadGroup&&) = delete;

    ~ThreadGroup() { stop(); }

    void stop() noexcept {
        {
            LockType l{m_mutex};
            m_stopped = true;
            m_condvar.notify_all();
        }
        std::for_each(m_threads.begin(), m_threads.end(),
                      [](std::thread& t) { t.join(); });
        m_threads.clear();
        m_stopped = false;
    }

    std::future<ResultType> post(TaskType&& t) {
        assert(t.valid());
        p_ensure_running();
        std::future<ResultType> res = t.get_future();
        {
            LockType l{m_mutex};
            m_tasks.emplace_back(std::move(t));
            m_condvar.notify_one();
        }
        return res;
    }

    template <typename TaskInputIterator, typename FutureOutputIterator>
    void post(TaskInputIterator task_begin, TaskInputIterator task_end,
              FutureOutputIterator out) {
        LockType l{m_mutex};
        p_ensure_running();
        for (; task_begin != task_end; ++task_begin) {
            std::future<ResultType> res = task_begin->get_future();
            m_tasks.emplace_back(std::move(*task_begin));
            *out = std::move(res);
            ++out;
        }
        m_condvar.notify_all();
    }

    template<typename RandomAccessIterator, typename Callable>
    void parallel_foreach(RandomAccessIterator begin, RandomAccessIterator end, Callable&& callable)
    {
        parallel_foreach_iterator(begin, end, [&] (const RandomAccessIterator& iter) {
            std::forward<Callable>(callable)(*iter);
        });
    }

    template<typename RandomAccessIterator, typename Callable>
    void parallel_foreach_iterator(RandomAccessIterator begin, RandomAccessIterator end, Callable&& callable)
    {
        if(begin == end) return;
        RandomAccessIterator last_beg = begin;
        std::vector<std::future<ResultType>> &futures = p_make_futures();
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(begin, end, num_threads() + 1,
                [&] (RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if(s_beg == end) return;
                    if(s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back([s_beg,s_end,&callable] (ThreadIndex) {
                            for(auto s_cur = s_beg; s_cur != s_end; ++s_cur) {
                                callable(s_cur);
                            }
                        });
                        futures.emplace_back(std::move(m_tasks.back().get_future()));
                    }
                }
            );
            m_condvar.notify_all();
        }
        for(RandomAccessIterator iter = last_beg; iter != end; ++iter) {
            callable(iter);
        }
        std::for_each(futures.begin(), futures.end(), [] (std::future<ResultType>& f) { f.get(); });
		futures.clear();
    }

    template<typename RandomAccessIterator, typename ContextType, typename Callable>
    void parallel_foreach(RandomAccessIterator begin, RandomAccessIterator end, const ContextType& context, Callable&& callable)
    {
        parallel_foreach_iterator(begin, end, context, [&] (ContextType& local_ctx, RandomAccessIterator iter) {
            std::forward<Callable>(callable)(local_ctx, *iter);
        });
    }

    template<typename RandomAccessIterator, typename ContextType, typename Callable>
    void parallel_foreach_iterator(RandomAccessIterator begin, RandomAccessIterator end, const ContextType& context, Callable&& callable)
    {
        if(begin == end) return;
        std::vector<std::future<ResultType>> &futures = p_make_futures();
        RandomAccessIterator last_beg = begin;
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(begin, end, num_threads() + 1,
                [&] (RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if(s_beg == end) return;
                    if(s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back([s_beg,s_end,&callable,&context] (ThreadIndex) {
                            ContextType local_context{context};
                            for(auto s_cur = s_beg; s_cur != s_end; ++s_cur) {
                                callable(local_context, s_cur);
                            }
                        });
                        futures.emplace_back(std::move(m_tasks.back().get_future()));
                    }
                }
            );
            m_condvar.notify_all();
        }
        ContextType local_context{context};
        for(RandomAccessIterator iter = last_beg; iter != end; ++iter) {
            callable(local_context, iter);
        }
        std::for_each(futures.begin(), futures.end(), [] (std::future<ResultType>& f) { f.get(); });
		futures.clear();
    }

    template<typename RandomAccessIterator,
             typename ContextFunction/*(ThreadIndex,RandomAccessIterator)*/,
             typename Callable,/*(Context&, ThreadIndex, RandomAccessIterator)*/
             typename FinalContextFunction/*(Context&, ThreadIndex)*/>
    void context_function_parallel_foreach_iterator(RandomAccessIterator begin, RandomAccessIterator end,
                                                    ContextFunction&& context_fn, Callable&& callable,
                                                    FinalContextFunction&& final_ctx_fn)
    {
        if(begin == end) return;
        std::vector<std::future<ResultType>> &futures = p_make_futures();
        RandomAccessIterator last_beg = begin;
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(begin, end, num_threads() + 1,
                [&] (RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if(s_beg == end) return;
                    if(s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back(
                            [s_beg,s_end,&context_fn,&callable,&final_ctx_fn] (ThreadIndex index) {
                                auto s_cur = s_beg;
                                auto& local_context = context_fn(index, s_cur);
                                for(++s_cur; s_cur != s_end; ++s_cur) {
                                    callable(local_context, index, s_cur);
                                }
                                final_ctx_fn(local_context, index);
                            }
                        );
                        futures.emplace_back(std::move(m_tasks.back().get_future()));
                    }
                }
            );
            m_condvar.notify_all();
        }
        auto& local_context = context_fn(0, last_beg);
        for(++last_beg; last_beg != end; ++last_beg) {
            callable(local_context, 0, last_beg);
        }
        final_ctx_fn(local_context, 0);
        std::for_each(futures.begin(), futures.end(), [] (std::future<ResultType>& f) { f.get(); });
		futures.clear();
    }

    template<typename Callable>
    void run_n_copies(std::size_t n, Callable&& callable) {
		std::vector<std::future<ResultType>> &futures = p_make_futures();
        futures.reserve(n - 1);
        {
            LockType l{m_mutex};
            p_ensure_running();
            for(std::size_t i = 0; i < n - 1; ++i) {
                m_tasks.emplace_back([&callable] (ThreadIndex) { std::forward<Callable>(callable)(); });
                futures.emplace_back(std::move(m_tasks.back().get_future()));
            }
            m_condvar.notify_all();
        }
		std::forward<Callable>(callable)();
		std::for_each(futures.begin(), futures.end(), [] (std::future<ResultType>& f) { f.get(); });
		futures.clear();
    }

    std::size_t num_threads() const noexcept { return m_thread_count; }

    std::size_t num_threads_running() const noexcept {
        return m_threads.size();
    }

  private:
    using LockType = std::unique_lock<std::mutex>;

    std::vector<std::future<ResultType>> &p_make_futures() {
        m_futures_cache.reserve(num_threads());
        return m_futures_cache;
    }

    void p_worker_main() {
		const ThreadIndex my_index = m_thread_index_counter++;
        for (;;) {
            TaskType task;
            {
                LockType l{m_mutex};
                if (!p_await_task(l, task))
                    return;
            }
            task(my_index);
        }
    }

    void p_ensure_running() {
        while (m_threads.size() < m_thread_count) {
            m_threads.emplace_back([this]() { p_worker_main(); });
        }
    }

    bool p_should_wake() const {
        if (m_stopped)
            return true;
        return !m_tasks.empty();
    }

    bool p_await_task(LockType& held_lock, TaskType& out) {
        m_condvar.wait(held_lock, [this]() -> bool { return p_should_wake(); });
        if (m_tasks.empty())
            return false;
        out = std::move(m_tasks.back());
        m_tasks.pop_back();
        return true;
    }

    std::vector<std::thread> m_threads;
    std::size_t m_thread_count;
    std::mutex m_mutex;
    std::condition_variable m_condvar;
    std::vector<TaskType> m_tasks;
	std::vector<std::future<ResultType>> m_futures_cache;
    bool m_stopped = false;
	std::atomic<ThreadIndex> m_thread_index_counter{1};
};

/**
 * A buffer used by ResultToLocalContext to store
 * each thread's local context objects in.
 */
template<typename ResultType>
class ResultToLocalContextBuffer
{
  public:
    template<typename... ResultConsArgs>
    explicit ResultToLocalContextBuffer(ResultConsArgs&&... args)
    {
        m_results.emplace_back(std::forward<ResultConsArgs>(args)...);
    }

    ResultType& front() noexcept { return m_results.front(); }
    const ResultType& front() const noexcept { return m_results.front(); }

    ResultType* create_pointer() {
        std::unique_lock<std::mutex> l{m_mutex};
        m_results.emplace_back(m_results.front());
        return &m_results.back();
    }

    template<typename ResultInitType, typename ReduceFn>
    ResultType reduce(ResultInitType&& init, ReduceFn&& reduce) {
        std::unique_lock<std::mutex> l{m_mutex};
        return std::reduce(std::next(m_results.begin()), m_results.end(),
                           std::forward<ResultInitType>(init), std::forward<ReduceFn>(reduce));
    }

    template<typename ForEachFn> void for_each(ForEachFn&& callback) {
        std::for_each(std::next(m_results.begin()), m_results.end(), std::forward<ForEachFn>(callback));
    }

  private:
    std::mutex m_mutex;
    std::list<ResultType> m_results;
};

/**
 * A local context structure that creates a new copy
 * of a given object in a std::list for each copy
 * made of the local context structure itself.
 * This is useful to allow a void ThreadGroup to
 * return values that need to be joined/reduced to
 * a single result later on.
 */
template<typename ResultType>
class ResultToLocalContext 
{
  public:
    explicit ResultToLocalContext(ResultToLocalContextBuffer<ResultType>& buffer) noexcept :
        m_buffer(&buffer),
        m_local_pointer(&buffer.front())
    {}

    explicit ResultToLocalContext(const ResultToLocalContext<ResultType>& other) noexcept :
        m_buffer(other.m_buffer),
        m_local_pointer(m_buffer->create_pointer())
    {}

    ResultType &get() noexcept { return *m_local_pointer; }
    const ResultType& get() const noexcept { return *m_local_pointer; }

  private:
    ResultToLocalContextBuffer<ResultType>* m_buffer;
    ResultType* m_local_pointer;
};

}

#endif
