#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <iostream>
#include <sammy/process.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace sammy {

class Process::Impl {
  public:
    Impl(const boost::filesystem::path& path,
         const std::vector<std::string>& args) {
        std::vector<const char*> cargs;
        cargs.push_back(path.c_str());
        for (const auto& arg : args) {
            cargs.push_back(arg.c_str());
        }
        cargs.push_back(nullptr);
        int proc_in_pipe[2];
        int proc_out_pipe[2];
        int e = ::pipe(&proc_in_pipe[0]);
        if (e < 0) {
            throw std::runtime_error("Failed to create process stdin pipe.");
        }
        e = ::pipe(&proc_out_pipe[0]);
        if (e < 0) {
            ::close(proc_in_pipe[0]);
            ::close(proc_in_pipe[1]);
            throw std::runtime_error("Failed to create process stdout pipe.");
        }
        pid_t pid = ::fork();
        if (pid < 0) {
            ::close(proc_in_pipe[0]);
            ::close(proc_in_pipe[1]);
            ::close(proc_out_pipe[0]);
            ::close(proc_out_pipe[1]);
            throw std::runtime_error("Failed to fork process.");
        }
        if (pid == 0) {
            // child process
            ::close(proc_in_pipe[1]);
            ::close(proc_out_pipe[0]);
            if (::dup2(proc_in_pipe[0], STDIN_FILENO) < 0) {
                std::abort();
            }
            if (::dup2(proc_out_pipe[1], STDOUT_FILENO) < 0) {
                std::abort();
            }
            ::close(proc_in_pipe[0]);
            ::close(proc_out_pipe[1]);
            ::execvp(cargs[0], const_cast<char* const*>(cargs.data()));
            std::abort();
        } else {
            // parent process
            m_pid = pid;
            m_proc_in_pipe = proc_in_pipe[1];
            m_proc_out_pipe = proc_out_pipe[0];
            ::close(proc_in_pipe[0]);
            ::close(proc_out_pipe[1]);
        }
    }

    ~Impl() {
        if (m_pid >= 0) {
            if (!m_exit_code) {
                if (!m_killed) {
                    terminate();
                }
                std::cerr << "Process not waited on before destruction."
                          << std::endl;
                std::abort();
            }
        }
        if (m_proc_in_pipe >= 0) {
            ::close(m_proc_in_pipe);
            m_proc_in_pipe = -1;
        }
        if (m_proc_out_pipe >= 0) {
            ::close(m_proc_out_pipe);
            m_proc_out_pipe = -1;
        }
    }

    void write_to_process(const char* data, std::size_t size) {
        std::size_t total_written = 0;
        while (total_written < size) {
            auto ssize = ::write(m_proc_in_pipe,
                                 static_cast<const void*>(data + total_written),
                                 size - total_written);
            if (ssize < 0) {
                throw std::runtime_error("Failed to write to process.");
            }
            total_written += static_cast<std::size_t>(ssize);
        }
    }

    void send_eof() {
        if (m_proc_in_pipe >= 0) {
            ::close(m_proc_in_pipe);
            m_proc_in_pipe = -1;
        }
    }

    std::size_t read_from_process(char* data, std::size_t size) {
        std::size_t total_read = 0;
        while (total_read < size) {
            auto ssize =
                ::read(m_proc_out_pipe, static_cast<void*>(data + total_read),
                       size - total_read);
            if (ssize < 0) {
                throw std::runtime_error("Failed to read from process.");
            }
            if (ssize == 0) {
                // eof reached
                return total_read;
            }
            total_read += static_cast<std::size_t>(ssize);
        }
        return total_read;
    }

    int exit_code() const { return m_exit_code.value(); }

    bool have_waited() const noexcept { return m_exit_code.has_value(); }

    void terminate() {
        if (m_pid != -1) {
            if (::kill(m_pid, SIGKILL) < 0) {
                throw std::runtime_error("Failed to terminate process.");
            }
            m_killed = true;
        }
    }

    int wait() {
        if (m_pid < 0) {
            throw std::runtime_error("Process not started.");
        }
        if (m_exit_code) {
            return *m_exit_code;
        }
        for (;;) {
            int status;
            int result = ::waitpid(m_pid, &status, 0);
            if (result < 0) {
                if (errno == EINTR) {
                    continue;
                }
                throw std::runtime_error("Failed to wait on process.");
            }
            if (WIFEXITED(status)) {
                m_exit_code = WEXITSTATUS(status);
                return *m_exit_code;
            }
            if (WIFSIGNALED(status)) {
                m_exit_code = -WTERMSIG(status);
                return *m_exit_code;
            }
        }
    }

  private:
    using ExitCode = std::optional<int>;

    pid_t m_pid{-1};         //< the process ID
    int m_proc_in_pipe{-1};  //< the write-end of the process stdin
    int m_proc_out_pipe{-1}; //< the read-end of the process stdout
    bool m_killed{false};    //< true if we have killed the process
    ExitCode m_exit_code;    //< exit code if we have waited
};

Process::Process(const boost::filesystem::path& path,
                 const std::vector<std::string>& args)
    : m_impl(std::make_unique<Impl>(path, args)) {}

Process::~Process() = default;

/**
 * Write data to the process.
 * Throws an error on failure.
 */
void Process::write_to_process(const char* data, std::size_t size) {
    m_impl->write_to_process(data, size);
}

/**
 * Read data from the process.
 * Returns the number of bytes read;
 * only less than size on reaching EOF.
 */
std::size_t Process::read_from_process(char* data, std::size_t size) {
    return m_impl->read_from_process(data, size);
}

/**
 * Terminate the process, using
 * SIGKILL on unix-like systems and TerminateProcess on windows.
 */
void Process::terminate() { m_impl->terminate(); }

/**
 * Wait for the process to exit and grab its exit code.
 */
int Process::wait() { return m_impl->wait(); }

/**
 * Get the exit code of the process, if it has exited and
 * we have waited on it; otherwise, throws an error.
 */
int Process::exit_code() const { return m_impl->exit_code(); }

/**
 * Check if we have created and waited on the process.
 */
bool Process::have_waited() const noexcept {
    if (!m_impl) {
        return false;
    }
    return m_impl->have_waited();
}

void Process::send_eof() { m_impl->send_eof(); }

void Process::reset() { m_impl.reset(); }

Process::Process(Process&&) noexcept = default;
Process& Process::operator=(Process&&) noexcept = default;

Process::Process() noexcept : m_impl(nullptr) {}

} // namespace sammy
