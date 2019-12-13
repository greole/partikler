/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

#ifndef SPHDATASTRUCTURES_H
#define SPHDATASTRUCTURES_H

#include "cgal/CGALTYPEDEFS.hpp"
#include "cgal/CGALHelper.hpp"
#include "FileIO.hpp"
#include "Logger.hpp"
#include "Helper.hpp"
#include "SearchCubes.hpp"
#include "Object.hpp"
#include "Vec3.hpp"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <x86intrin.h>


// Move to NeighboursSearch
struct FirstLastPair{
  size_t first;
  size_t last;
};

struct VectorPair {
    // Stores the distance of particles on different
    // STL surfaces

    Vec3 on;
    Vec3 no;
};

struct NeighbourIdHalfStencil {
  // Stores the stride of a domain, ie the difference of search cubes ids
  // for the 26 neighbour cubes. Since the stride is symmetric
  // only half of the 26 values are stored. Addition of the stencil
  // values yields the upper and subtraction the lower neighbours
  // Values are Stored as size_t for simpler addition without type conversion.

  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |x  x  x      |1  2  3      |10 11  12 // back
  //  |x  x  x      |x  x  0      |7   8   9 // middle
  //  |x  x  x      |x  x  x      |4   5   6 // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3
  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |06 07 08     | x  x  x     | x  x  x // back
  //  |03 04 05     |12  x  x     | x  x  x // middle
  //  |00 01 02     |09 10 11     | x  x  x // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3

  std::vector<size_t> stencil;

  NeighbourIdHalfStencil (size_t nx, size_t ny) {
    // TODO leave it a size_t, iterate only till 12, since
    // the stencil is symmetric anyway
    const size_t ny_nx = nx*ny;

    stencil = std::vector<size_t> {
                                   1,              // right
                                   nx - 1,         // back left
                                   nx    ,         // back centre
                                   nx + 1,         // back right
                                   ny_nx - nx - 1, // upper front left
                                   ny_nx - nx    , // upper front centre
                                   ny_nx - nx + 1, // upper front right
                                   ny_nx       -1, // upper middle left
                                   ny_nx         , // upper middle centre
                                   ny_nx      + 1, // upper middle right
                                   ny_nx + nx - 1, // upper back left
                                   ny_nx + nx    , // upper back centre
                                   ny_nx + nx + 1  // upper back right
    };
  };
};


#endif
