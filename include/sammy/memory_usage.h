#ifndef BARRAGE_MEMORY_USAGE_H_INCLUDED_
#define BARRAGE_MEMORY_USAGE_H_INCLUDED_

#include <cstdint>
#include <cstddef>
#include <cstdlib>

namespace sammy {
    
/**
 * @brief Get the current peak resident set size (RSS) in bytes,
 * i.e., the maximum RSS since the start of the process, so far.
 */
inline std::size_t current_peak_rss();

}

#if defined(__APPLE__) && defined(__MACH__)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if(getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss;
}

#elif __has_include(<sys/resource.h>)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if(getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss * 1024;
}

#elif __has_include(<psapi.h>) && __has_include(<windows.h>) && defined(WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>

std::size_t sammy::current_peak_rss() {
    PROCESS_MEMORY_COUNTERS pmc;
    if(!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        std::abort();
    }
    return pmc.PeakWorkingSetSize;
}

#endif

#endif
