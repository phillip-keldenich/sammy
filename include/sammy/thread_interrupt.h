#ifndef SAMMY_THREAD_INTERRUPT_H_INCLUDED_
#define SAMMY_THREAD_INTERRUPT_H_INCLUDED_

#include <atomic>
#include <exception>
#include <stdexcept>

namespace sammy {

class InterruptError : public std::exception {
  public:
    const char* what() const noexcept override { return "Interrupted"; }
};

inline std::atomic<bool>*& get_interrupt_flag_ptr() noexcept {
    thread_local std::atomic<bool>* flag = nullptr;
    return flag;
}

inline void set_interrupt_flag_ptr(std::atomic<bool>* flag) noexcept {
    get_interrupt_flag_ptr() = flag;
}

inline bool peek_interrupt_flag() noexcept {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        return false;
    return iflag->load();
}

inline bool get_and_clear_interrupt_flag() noexcept {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        return false;
    bool fv = iflag->load();
    if (!fv)
        return false;
    fv = iflag->exchange(false);
    return fv;
}

inline void trigger_interrupt() {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        throw std::logic_error("Interrupt flag not set!");
    iflag->store(true);
}

inline void throw_if_interrupted() {
    if (get_and_clear_interrupt_flag()) {
        throw InterruptError{};
    }
}

} // namespace sammy

#endif
