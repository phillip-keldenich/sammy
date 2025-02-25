#ifndef SAMMY_BASIC_CUDA_ITERATION_H_INCLUDED_
#define SAMMY_BASIC_CUDA_ITERATION_H_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>

/**
 * This file exists as a workaround to bugs with nvcc (or gcc/clang/stdlibc++,
 * depends on who you ask) which causes compilation problems when nvcc
 * is exposed to some standard C++ headers.
 */

/**
 * Conditional compilation (if no CUDA support,
 * the methods in this file cannot be used and
 * will lead to linker errors, but the file can be included).
 */
#ifndef SAMMY_CUDA_SUPPORTED
using cudaError_t = int;
#define SAMMY_HD
#else
#include <cuda_runtime.h>
#define SAMMY_HD __host__ __device__
#endif
#include <atomic>


namespace sammy {

/**
 * Exception thrown when CUDA errors occur.
 */
class CUDAError : public std::runtime_error {
  public:
    explicit CUDAError(cudaError_t err);
    explicit CUDAError(cudaError_t err, const std::string& info);

    /**
     * Get the CUDA error code.
     */
    cudaError_t error_code() const noexcept;

    /**
     * Get the CUDA error string.
     */
    const char* error_string() const noexcept;

    /**
     * Get the CUDA error name.
     */
    const char* error_name() const noexcept;

  private:
    cudaError_t m_error_code;
};

/**
 * CUDA use mode: if available, use CUDA; on error, fallback to CPU.
 */
static constexpr int CUDA_USAGE_IF_AVAILABLE = 0;

/**
 * CUDA use mode: use CUDA; on error, fail.
 */
static constexpr int CUDA_USAGE_FORCED = 1;

/**
 * CUDA use mode: do not use CUDA.
 */
static constexpr int CUDA_USAGE_DISABLED = 2;

/**
 * CUDA use mode: had an error previously, don't use CUDA.
 */
static constexpr int CUDA_USAGE_HAD_ERROR = 3;


#ifdef SAMMY_CUDA_SUPPORTED

/**
 * Usage info flag for CUDA.
 */
inline std::atomic<int>& cuda_usage_info() {
    static std::atomic<int> usage_info{CUDA_USAGE_IF_AVAILABLE};
    return usage_info;
}

inline bool should_use_cuda() {
    switch(cuda_usage_info()) {
        default:
        case CUDA_USAGE_IF_AVAILABLE:
        case CUDA_USAGE_FORCED:
            return true;

        case CUDA_USAGE_HAD_ERROR:
        case CUDA_USAGE_DISABLED:
            return false;
    }
}

inline void had_cuda_error_(const CUDAError& error) {
    auto& info = cuda_usage_info();
    int loaded = info.load();
    if(loaded == CUDA_USAGE_FORCED) {
        throw error;
    }
    info.store(CUDA_USAGE_HAD_ERROR);
}

inline void had_cuda_error(const CUDAError& error) noexcept {
    had_cuda_error_(error);
}

inline void set_cuda_mode(int mode) {
    cuda_usage_info().store(mode);
}

#else

inline constexpr bool should_use_cuda() { return false; }

inline void had_cuda_error(const CUDAError& error) noexcept {}

inline void set_cuda_mode(int mode) {
    if(mode == CUDA_USAGE_FORCED) {
        throw std::runtime_error("Requested CUDA mode 'force' without CUDA support compiled in!");
    }
}

#endif

namespace detail {

void cuda_malloc_or_throw(void** ptr, std::size_t bytes);
void cuda_memcpy_htd_or_throw(void* dst, const void* src, std::size_t bytes);
void cuda_memcpy_dth_or_throw(void* dst, const void* src, std::size_t bytes);
void cuda_free(void* ptr);

/**
 * We produce rows of output.
 * Each row is a bitset of size num_literals.
 * Each bitset consists of 32-bit integers.
 * THREADS_PER_ROW threads produce one output bitset together.
 */
constexpr static SAMMY_HD std::size_t BLOCK_SIZE() { return 256; }
constexpr static SAMMY_HD std::size_t THREADS_PER_ROW() { return 16; }
constexpr static SAMMY_HD std::size_t ROWS_PER_BLOCK() { return BLOCK_SIZE() / THREADS_PER_ROW(); }
constexpr static SAMMY_HD std::size_t GOAL_ROWS_PER_CALL() { return 256 * ROWS_PER_BLOCK(); }

void cuda_call_bit_filter_kernel(const std::uint32_t* bit_data, std::size_t u32_per_bitset,
                                 const std::uint32_t* classes_with_literal, std::size_t num_literals,
                                 const std::uint32_t* classes_with_literal_offsets,
                                 std::uint32_t* output_buffer, std::size_t begin_row, std::size_t num_rows);

void call_cuda_extract_kernel(const std::uint32_t* device_bit_data, std::size_t u32_per_bitset,
                              std::uint32_t* output_buffer, std::size_t num_literals,
                              std::size_t begin_var, std::size_t num_vars, std::size_t num_classes);

}

}

#endif

