#include <sammy/basic_cuda_iteration.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

namespace sammy {

CUDAError::CUDAError(cudaError_t err) :
    std::runtime_error(
        std::string("CUDA error ") + 
        cudaGetErrorName(err) + ": " +
        cudaGetErrorString(err)),
    m_error_code(err)
{}

CUDAError::CUDAError(cudaError_t err, const std::string& info) :
    std::runtime_error(
        "CUDA error (" + info + ") " +
        cudaGetErrorName(err) + ": " +
        cudaGetErrorString(err)),
    m_error_code(err)
{}

cudaError_t CUDAError::error_code() const noexcept {
    return m_error_code;
}

const char* CUDAError::error_string() const noexcept {
    return cudaGetErrorString(m_error_code);
}

const char* CUDAError::error_name() const noexcept {
    return cudaGetErrorName(m_error_code);
}


namespace detail {

void cuda_malloc_or_throw(void** ptr, std::size_t bytes) {
    cudaError_t err = cudaMalloc(ptr, bytes);
    if(err != cudaSuccess) {
        throw CUDAError(err, "cudaMalloc(" + std::to_string(bytes) + " bytes)");
    }
}

void cuda_memcpy_htd_or_throw(void* dst, const void* src, std::size_t bytes) {
    cudaError_t err = cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
        throw CUDAError(err, "cudaMemcpy(" + std::to_string(bytes) + " bytes, host -> device)");
    }
}

void cuda_memcpy_dth_or_throw(void* dst, const void* src, std::size_t bytes) {
    cudaError_t err = cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
        throw CUDAError(err, "cudaMemcpy(" + std::to_string(bytes) + " bytes, device -> host)");
    }
}

void cuda_free(void* ptr) {
    cudaFree(ptr);
}

static void __global__ 
    cuda_bit_filter_kernel(const std::uint32_t* bit_data, std::size_t u32_per_bitset,
                           const std::uint32_t* classes_with_literal, std::size_t num_literals,
                           const std::uint32_t* classes_with_literal_offsets,
                           std::uint32_t* output_buffer, std::size_t begin_row, std::size_t num_rows)
{
    std::size_t block_idx(blockIdx.x);
    std::size_t thread_idx(threadIdx.x);
    std::size_t row_offs = block_idx * ROWS_PER_BLOCK() + thread_idx / THREADS_PER_ROW();
    std::uint32_t intra_row_offs = thread_idx % THREADS_PER_ROW();
    if(row_offs >= num_rows) { return; }
    std::size_t row_idx = begin_row + row_offs;
    if(row_idx >= num_literals) { 
        __trap(); // cause abnormal termination if we were invoked incorrectly
    }
    const std::uint32_t* cwl_beg = &classes_with_literal[classes_with_literal_offsets[row_idx]];
    const std::uint32_t* cwl_end = &classes_with_literal[classes_with_literal_offsets[row_idx + 1]];
    std::uint32_t* output_row = &output_buffer[row_offs * u32_per_bitset];
    for(std::size_t init_idx = intra_row_offs; init_idx < u32_per_bitset; init_idx += THREADS_PER_ROW()) {
        output_row[init_idx] = 0;
    }
    for(; cwl_beg != cwl_end; ++cwl_beg) {
        std::uint32_t current_class = *cwl_beg;
        const std::uint32_t* bitset_beg = &bit_data[current_class * u32_per_bitset];
        for(std::size_t idx = row_idx / 32 + intra_row_offs; idx < u32_per_bitset; idx += THREADS_PER_ROW()) {
            output_row[idx] |= bitset_beg[idx];
        }
    }
}

void cuda_call_bit_filter_kernel(const std::uint32_t* bit_data, std::size_t u32_per_bitset,
                                 const std::uint32_t* classes_with_literal, std::size_t num_literals,
                                 const std::uint32_t* classes_with_literal_offsets,
                                 std::uint32_t* output_buffer, std::size_t begin_row, std::size_t num_rows)
{
    auto num_blocks = num_rows / ROWS_PER_BLOCK();
	if(num_rows % ROWS_PER_BLOCK()) ++num_blocks;
    cuda_bit_filter_kernel<<<num_blocks, BLOCK_SIZE()>>>(
        bit_data, u32_per_bitset,
        classes_with_literal, num_literals,
        classes_with_literal_offsets,
        output_buffer, begin_row, num_rows);
	auto err = cudaPeekAtLastError();
	if(err != cudaSuccess) {
		throw CUDAError(err, "CUDA kernel launch error (bit filter)");
	}
	err = cudaDeviceSynchronize();
	if(err != cudaSuccess) {
		throw CUDAError(err, "CUDA kernel asynchronous error (bit filter)");
	}
}

static void __global__ cuda_extract_kernel(const std::uint32_t* device_bit_data, std::size_t u32_per_bitset,
                                           std::uint32_t* output_buffer, std::size_t num_literals,
                                           std::size_t begin_var, std::size_t num_vars, std::size_t num_classes)
{
    std::size_t block_idx(blockIdx.x);
    std::size_t thread_idx(threadIdx.x);
    std::size_t var_offs = block_idx * ROWS_PER_BLOCK() + thread_idx / THREADS_PER_ROW();
    std::uint32_t intra_row_offs = thread_idx % THREADS_PER_ROW();
    if(var_offs >= num_vars) { return; }
    std::size_t var_idx = begin_var + var_offs;
    std::size_t pos_literal_index = 2 * var_idx;
    std::size_t bitset_wrd = pos_literal_index / 32;
    std::size_t bitset_bit = pos_literal_index % 32;
    if(pos_literal_index + 1 >= num_literals) {
        __trap(); // cause abnormal termination if we were invoked incorrectly
    }
    std::uint32_t* pos_out_row = &output_buffer[2 * var_offs * u32_per_bitset];
    std::uint32_t* neg_out_row = &output_buffer[(2 * var_offs + 1) * u32_per_bitset];
    for(std::size_t init_idx = intra_row_offs; init_idx < u32_per_bitset; init_idx += THREADS_PER_ROW()) {
        pos_out_row[init_idx] = 0;
        neg_out_row[init_idx] = 0;
    }
    for(std::size_t ci = 0; ci < num_classes; ++ci) {
        const std::uint32_t* class_begin = &device_bit_data[ci * u32_per_bitset];
        std::uint32_t is_positive = class_begin[bitset_wrd] & (1 << bitset_bit);
        std::uint32_t* out_ptr = is_positive ? pos_out_row : neg_out_row;
        for(std::size_t i = intra_row_offs; i < u32_per_bitset; i += THREADS_PER_ROW()) {
            out_ptr[i] |= class_begin[i];
        }
    }
}

void call_cuda_extract_kernel(const std::uint32_t* device_bit_data, std::size_t u32_per_bitset,
                              std::uint32_t* output_buffer, std::size_t num_literals,
                              std::size_t begin_var, std::size_t num_vars, std::size_t num_classes)
{
    auto num_blocks = num_vars / ROWS_PER_BLOCK();
    if(num_vars % ROWS_PER_BLOCK()) ++num_blocks;
    cuda_extract_kernel<<<num_blocks, BLOCK_SIZE()>>>(
        device_bit_data, u32_per_bitset,
        output_buffer, num_literals,
        begin_var, num_vars, num_classes);
    auto err = cudaPeekAtLastError();
    if(err != cudaSuccess) {
        throw CUDAError(err, "CUDA kernel launch error (extract)");
    }
    err = cudaDeviceSynchronize();
    if(err != cudaSuccess) {
        throw CUDAError(err, "CUDA kernel asynchronous error (extract)");
    }
}

}
}
