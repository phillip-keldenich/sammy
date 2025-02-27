#ifndef SAMMY_ERROR_H_INCLUDED_
#define SAMMY_ERROR_H_INCLUDED_

#include <exception>
#include <stdexcept>
#include <string>

namespace sammy {

class UNSATError : public std::runtime_error {
  public:
    UNSATError() : std::runtime_error("The given model is unsatisfiable") {}
};

} // namespace sammy

#endif
