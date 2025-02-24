# External SAT solvers
This directory contains the raw source of external SAT solvers
used by our algorithms. They are included in several ways, depending
on their build system; generally, we integrate the building of external
SAT solvers into our build instead of relying on pre-compiled binaries
or bringing them in as external dependencies.

Each of the solvers is then wrapped in a binding that makes the solvers
accessible in a uniform way, allowing the use of template parameters
to exchange one external SAT solver for another.

