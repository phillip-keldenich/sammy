# Conanfile for Gurobi

This is a Conan package for the Gurobi Optimizer. Gurobi is a commercial optimization solver, and you need to have a license to use it. You can get a free academic license from the [Gurobi website](https://www.gurobi.com/).

You have to install this package manually as it is not available in the Conan Center Index (because it is a commercial software, and we are not its owner/maintainer).
This package was developed by us to help you to create your own packages that depend on Gurobi, but we have no affiliation with Gurobi.

You can make gurobi available to your environment by running

```bash
conan create . -pr default -s build_type=Release --build=missing
```

Afterwards, you can use it in your own projects by adding `gurobi/11.0.0` to your `requires` field in your `conanfile.txt` or `conanfile.py`.

```
[requires]
gurobi/[>=11.0.0]

[generators]
CMakeDeps
CMakeToolchain
```