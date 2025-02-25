#ifndef SAMMY_ERROR_H_INCLUDED_
#define SAMMY_ERROR_H_INCLUDED_

#include <string>
#include <stdexcept>
#include <exception>

namespace sammy {

class UNSATError : public std::runtime_error {
public:
    UNSATError() : std::runtime_error("The given model is unsatisfiable") {}
};

}

#endif
