#ifndef SAMMY_PROCESS_H_INCLUDED_
#define SAMMY_PROCESS_H_INCLUDED_

#include <boost/filesystem.hpp>
#include <memory>

namespace sammy {

/**
 * A very simplistic process management class,
 * since Boost.Process is a complete mess when it comes to managing and
 * waiting on processes.
 * It can't even wait on a process in a thread-safe manner,
 * even without any time limits, and with completely unrelated
 * boost::process::child objects.
 */
class Process {
  public:
    Process() noexcept;

    Process(const boost::filesystem::path& path,
            const std::vector<std::string>& args);

    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;

    Process(Process&&) noexcept;
    Process& operator=(Process&&) noexcept;

    /**
     * Destructor aborts if the process
     * is still not waited for at destruction time.
     */
    ~Process();

    /**
     * Write data to the process.
     * Throws an error on failure.
     */
    void write_to_process(const char* data, std::size_t size);

    /**
     * Send EOF to the process.
     */
    void send_eof();

    /**
     * Read data from the process.
     * Returns the number of bytes read;
     * only less than size on reaching EOF.
     */
    std::size_t read_from_process(char* data, std::size_t size);

    /**
     * Terminate the process, using
     * SIGKILL on unix-like systems and TerminateProcess on windows.
     */
    void terminate();

    /**
     * Wait for the process to exit and grab its exit code.
     */
    int wait();

    /**
     * Get the exit code of the process, if it has exited and
     * we have waited on it; otherwise, throws an error.
     */
    int exit_code() const;

    /**
     * Check if the process has ever been started.
     */
    bool valid() const noexcept { return m_impl != nullptr; }

    /**
     * Reset the process to represent not-a-process.
     */
    void reset();

    /**
     * Check if we have created and waited on the process.
     */
    bool have_waited() const noexcept;

  private:
    class Impl;

    std::unique_ptr<Impl> m_impl;
};

} // namespace sammy

#endif
