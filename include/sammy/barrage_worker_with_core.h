#ifndef SAMMY_BARRAGE_WORKER_WITH_CORE_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_WITH_CORE_H_INCLUDED_

#include "barrage.h"
#include "thread_interrupt.h"

namespace sammy {

template <typename CoreType>
class PortfolioElementWithCore : public PortfolioElement {
  public:
    using CoreFactoryType = std::function<std::unique_ptr<CoreType>(
        PortfolioSolver*, PortfolioElementWithCore*)>;

    PortfolioElementWithCore(PortfolioSolver* solver,
                             CoreFactoryType core_factory,
                             const std::string& description)
        : PortfolioElement(solver), m_core_factory(std::move(core_factory)),
          m_description(description) {}

    EventRecorder* get_mutable_recorder() { return &m_recorder; }

    std::mutex& get_mutex() noexcept { return mutex; }

  protected:
    void main() override {
        if (should_terminate.load()) {
            return;
        }
        set_interrupt_flag_ptr(&should_terminate);
        try {
            p_construct_core();
        } catch (const InterruptError&) {
            return;
        }
        if (should_terminate.load()) {
            return;
        }
        m_core->main();
    }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) override {
        if (m_core) {
            if (should_terminate.load()) {
                m_core->termination_flag_set();
            } else {
                m_core->interrupt_if_necessary(info);
            }
        }
        // otherwise, we are constructing the core and,
        // if termination is wanted, should_terminate is automatically set
        // which causes the thread to interrupt
        events = 0;
    }

    const EventRecorder* get_recorder() const override { return &m_recorder; }

    /**
     * If this element has any, synchronize the event recorder.
     */
    virtual void synchronize_recorder(const EventRecorder& other) override {
        m_recorder.synchronize_with(other);
    }

    /**
     * If this element has any, set the quiet flag of the event recorder.
     */
    virtual void set_recorder_quiet(bool quiet) override {
        m_recorder.set_print_events(!quiet);
    }

    /**
     * Get a description of this element (it it has any).
     */
    virtual std::string get_description() const override {
        return m_description;
    }

    std::unique_ptr<CoreType> m_core;
    CoreFactoryType m_core_factory;
    EventRecorder m_recorder;
    std::string m_description;

  private:
    void p_construct_core() {
        std::unique_ptr<CoreType> core = m_core_factory(this->solver, this);
        {
            std::unique_lock l{mutex};
            m_core = std::move(core);
        }
    }
};

} // namespace sammy

#endif
