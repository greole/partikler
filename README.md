![](https://api.travis-ci.org/greole/partikler.svg?branch=master)

# Build Instructions

    mkdir build && cd build
    cmake ..
    make

## Dependencies

The partikler depends on CGAL and google tests. You can use the build process to fetch the depencies for you using the following flags.

    cmake  -DBUILD_GOOGLE_TEST=TRUE -DBUILD_CGAL=TRUE ..



# Usage

    ./partikler input.stl output.off dx

