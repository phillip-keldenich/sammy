#ifndef SAMMY_LOGGING_H_INCLUDED_
#define SAMMY_LOGGING_H_INCLUDED_

#include "literals.h"

namespace sammy {

template <typename StreamType, typename ContainerType>
inline void print_clause_internal(StreamType& stream,
                                  const ContainerType& clause) {
    bool first = true;
    for (Lit l : clause) {
        if (!first)
            stream << ' ';
        stream << l;
        first = false;
    }
}

template <typename StreamType, typename ContainerType>
inline void print_clause_external(StreamType& stream,
                                  const ContainerType& clause) {
    bool first = true;
    for (Lit l : clause) {
        if (!first)
            stream << ' ';
        stream << lit::externalize(l);
        first = false;
    }
}

} // namespace sammy

#endif
