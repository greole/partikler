<div align="center">
  <img src="https://raw.githubusercontent.com/greole/partikler/master/doc/LogoSmall.png"><br>
</div>

-----------------
[![Build Status](https://travis-ci.org/greole/partikler.svg?branch=master)](https://travis-ci.org/greole/partikler)
[![Coverage Status](https://coveralls.io/repos/github/greole/partikler/badge.svg?branch=master)](https://coveralls.io/github/greole/partikler?branch=master)
[![Documentation](https://codedocs.xyz/greole/partikler.svg)](https://codedocs.xyz/greole/partikler/)

# Features

 - Generate SPH particle boundary from .stl files
 - Expression templates for simple but efficient transport equation implementation

<div align="center">
  <img src="https://raw.githubusercontent.com/greole/partikler/master/doc/Screenshot.png"><br>
</div>

# Models

## Kernel
 - Gaussian
 - Wendland

## Surface Tension
 - Akinci





# Build Instructions

    mkdir build && cd build
    cmake ..
    make

## Dependencies

The partikler depends on CGAL Version 3.13 and boost yap and google tests. You can use the build process to fetch the depencies for you using the following flags.

    cmake  -DBUILD_GOOGLE_TEST=TRUE -DBUILD_CGAL=TRUE -DBUILD_HDF5=ON ..

### Ubuntu

To install the depencies for build google tests and cgal you can use

    sudo apt install libboost-all-dev libyaml-cpp-dev libmpfr-dev  libgmp-dev libomp-dev


# Usage

    ./partikler --config=<path/to/config.yaml>

