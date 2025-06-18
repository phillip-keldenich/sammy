#ifndef BARRAGE_MEMORY_USAGE_H_INCLUDED_
#define BARRAGE_MEMORY_USAGE_H_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <cstdlib>

#if defined(__has_include) && defined(__linux__) && defined(__GNUC__)
#if __has_include(<malloc.h>)
#include <malloc.h>
#endif
#endif

namespace sammy {

/**
 * @brief Attempt to return memory allocated and freed by the process
 *        to the operating system.
 *
 * Currently, this function only does something on Linux with
 * `malloc_trim(0)` available, which can help reduce out-of-memory (OOM)
 * errors because malloc seems to very eagerly hold on to memory for
 * small allocations (which mostly occur during the initial phase),
 * resulting in huge amounts of memory being allocated uselessly,
 * triggering OOM kills after the initial phase.
 */
inline void return_memory_to_os() {
#if defined(__has_include) && defined(__linux__) && defined(__GNUC__)
#if __has_include(<malloc.h>)
    ::malloc_trim(0); // really helps with OOM errors on Linux
#endif
#endif
}

/**
 * @brief Get the current peak resident set size (RSS) in bytes,
 * i.e., the maximum RSS since the start of the process, so far.
 */
inline std::size_t current_peak_rss();

} // namespace sammy

#if defined(__APPLE__) && defined(__MACH__)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if (getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss;
}

#elif __has_include(<sys/resource.h>)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if (getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss * 1024;
}

#elif __has_include(<psapi.h>) && __has_include(<windows.h>) && defined(WIN32)

#define WIN32_LEAN_AND_MEAN
#include <psapi.h>
#include <windows.h>

std::size_t sammy::current_peak_rss() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        std::abort();
    }
    return pmc.PeakWorkingSetSize;
}

#endif

#endif
