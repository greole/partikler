<div align="center">
  <img src="https://raw.githubusercontent.com/greole/partikler/master/doc/LogoSmall.png"><br>
</div>

-----------------
![](https://api.travis-ci.org/greole/partikler.svg?branch=master)

# Build Instructions

    mkdir build && cd build
    cmake ..
    make

## Dependencies

The partikler depends on CGAL and google tests. You can use the build process to fetch the depencies for you using the following flags.

    cmake  -DBUILD_GOOGLE_TEST=TRUE -DBUILD_CGAL=TRUE ..

### Ubuntu

To install the depencies for build google tests and cgal you can use

    sudo apt install libboost-all-dev libyaml-cpp-dev libmpfr-dev  libgmp-dev libomp-dev


# Usage

    ./partikler --config=<path/to/config.yaml>

