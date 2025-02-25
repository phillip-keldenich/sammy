#ifndef SAMMY_CUDA_ITERATION_H_INCLUDED_
#define SAMMY_CUDA_ITERATION_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "basic_cuda_iteration.h"

namespace sammy {
namespace detail {

template<typename Container, 
         std::enable_if_t<sizeof(DynamicBitset::Block) == 0 * sizeof(Container) + 
                          sizeof(std::uint64_t), int> = 0>
static void to_prepare_buffer(DynamicBitset::Block b, Container& buffer) {
    buffer.push_back(std::uint32_t(b & std::numeric_limits<std::uint32_t>::max()));
    buffer.push_back(std::uint32_t(b >> 32));
}

template<typename Container, 
         std::enable_if_t<sizeof(DynamicBitset::Block) == 0 * sizeof(Container) + 
                          sizeof(std::uint32_t), int> = 0>
static void to_prepare_buffer(DynamicBitset::Block b, Container& buffer) {
    buffer.push_back(std::uint32_t(b));
}

template<typename T>
class CUDADevicePointer {
    static_assert(std::is_pod_v<T>, "T must be POD");

    std::remove_cv_t<T>* m_device_ptr{nullptr};
    std::size_t m_device_size{0};

  public:
    /**
     * Create an empty device pointer.
     */
    CUDADevicePointer() = default;

    /**
     * Create a device pointer pointing to a buffer
     * of size objects of type T.
     */
    CUDADevicePointer(std::size_t size) {
        std::size_t bytes = size * sizeof(T) * (CHAR_BIT / 8);
        cuda_malloc_or_throw((void**)(&m_device_ptr), bytes);
        m_device_size = size;
    }

    /**
     * Create a device pointer pointing to a buffer
     * of size objects of type T, and copy the data
     * from the given vector.
     */
    CUDADevicePointer(const std::vector<std::remove_cv_t<T>>& host_data) {
        std::size_t bytes = host_data.size() * sizeof(T) * (CHAR_BIT / 8);
        cuda_malloc_or_throw((void**)(&m_device_ptr), bytes);
        try {
            cuda_memcpy_htd_or_throw(static_cast<void*>(m_device_ptr), static_cast<const void*>(host_data.data()), bytes);
        } catch(...) {
            cuda_free(static_cast<void*>(m_device_ptr));
            m_device_ptr = nullptr;
            throw;
        }
        m_device_size = host_data.size();
    }

    /**
     * Destroy the object (freeing the device memory).
     */
    ~CUDADevicePointer() {
        if(m_device_ptr) {
            cuda_free(static_cast<void*>(m_device_ptr));
        }
    }

    /**
     * Get the pointer to the device memory.
     */
    T* get() const noexcept {
        return m_device_ptr;
    }

    /**
     * Copy the device memory to the given host buffer.
     */
    void copy_to_host(std::vector<std::remove_cv_t<T>>& host_buffer) const {
        host_buffer.resize(m_device_size);
        cuda_memcpy_dth_or_throw(host_buffer.data(), m_device_ptr, 
                                 m_device_size * sizeof(T) * (CHAR_BIT / 8));
    }

    /**
     * Copy the device memory to a vector and return it.
     */
    std::vector<std::remove_cv_t<T>> get_host_copy() const {
        std::vector<std::remove_cv_t<T>> result;
        copy_to_host(result);
        return result;
    }
};

}

}

#endif
