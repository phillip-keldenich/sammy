"""
This file instructs scitkit-build how to build the module. This is very close
to the classical setuptools, but wraps some things for us to build the native
modules easier and more reliable.

For a proper installation with `pip install .`, you additionally need a
`pyproject.toml` to specify the dependencies to load this `setup.py`.

You can use `python3 setup.py install` to build and install the package
locally, with verbose output. To build this package in place, which may be
useful for debugging, use `python3 setup.py develop`. This will build
the native modules and move them into your source folder.

The setup options are documented here:
https://scikit-build.readthedocs.io/en/latest/usage.html#setup-options
"""

from glob import glob

from setuptools import find_packages
from skbuild_conan import setup


def readme():
    # Simply return the README.md as string
    with open("README.md") as file:
        return file.read()


setup(  # https://scikit-build.readthedocs.io/en/latest/usage.html#setup-options
    # ~~~~~~~~~ BASIC INFORMATION ~~~~~~~~~~~
    name="sample-analyzer",
    description="",
    version="0.1.0",
    long_description=readme(),
    long_description_content_type="text/markdown",
    author="TU Braunschweig, IBR, Algorithms Group",
    author_email="TBD",
    classifiers=[
        "Development Status :: 4 - Beta",
        # "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    # ~~~~~~~~~~~~ CRITICAL PYTHON SETUP ~~~~~~~~~~~~~~~~~~~
    # This project structures defines the python packages in a subfolder.
    # Thus, we have to collect this subfolder and define it as root.
    packages=find_packages(
        "src", exclude=["tests"]
    ),  # Include all packages in `./python`.
    package_dir={"": "src"},  # The root for our python package is in `./src`.
    python_requires=">=3.8",  # lowest python version supported.
    data_files=[
    ],
    entry_points={"console_scripts": ["samplns=samplns.__main__:main"]},
    cmake_minimum_required_version="3.24",
)