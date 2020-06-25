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

#ifndef SEARCHCUBES_H
#define SEARCHCUBES_H

#include <omp.h>    // omp_get_num_threads
#include <stddef.h> // for size_t
#include <vector>   // for vector

#include <algorithm> // for std::sort

#include "Scalar.hpp"
#include "Vec3.hpp"

std::pair<Vec3, Vec3> bounding_box(std::vector<Vec3> const &particles);

struct SubDivision {
    // 3*4bytes = 12bytes
    // Max 65535**3 cubes, which are approx (65535*3)**3 particles
    unsigned int nx, ny, nz;
};

struct SearchCubeDomain {

    Vec3 min;
    Vec3 max;

    Scalar dx;
    Scalar idx;

    // 32 byte

    SubDivision n;
    size_t nt;

    // 24 byte

    /* size_t padding; */
};

template <class T>
size_t position_to_cube_id(SearchCubeDomain scd, const T &p) {
    // TODO test if speed up n_cubes are copied to a const size_t nx ...

    const size_t nx = (size_t)scd.n.nx;
    const size_t ny = (size_t)scd.n.ny;
    // const size_t nz = (size_t)scd.n.nz;

    // use idx instead of dx
    const Scalar idx = scd.idx;
    const Scalar mx = scd.min[0];
    const Scalar my = scd.min[1];
    const Scalar mz = scd.min[2];

    const size_t i {(size_t)((p[0] - mx) * idx)};
    const size_t j {(size_t)((p[1] - my) * idx)};
    const size_t k {(size_t)((p[2] - mz) * idx)};

    return i + nx * j + nx * ny * k;
}

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

    NeighbourIdHalfStencil(size_t nx, size_t ny) {
        // TODO leave it a size_t, iterate only till 12, since
        // the stencil is symmetric anyway
        const size_t ny_nx = nx * ny;

        stencil = std::vector<size_t> {
            1,              // right
            nx - 1,         // back left
            nx,             // back centre
            nx + 1,         // back right
            ny_nx - nx - 1, // upper front left
            ny_nx - nx,     // upper front centre
            ny_nx - nx + 1, // upper front right
            ny_nx - 1,      // upper middle left
            ny_nx,          // upper middle centre
            ny_nx + 1,      // upper middle right
            ny_nx + nx - 1, // upper back left
            ny_nx + nx,     // upper back centre
            ny_nx + nx + 1  // upper back right
        };
    };
};

SearchCubeDomain
initSearchCubeDomain(const std::pair<Vec3, Vec3> &bound_box, Scalar dx);

struct NeighbourPair {
    size_t ownId;
    size_t neighId;
};

struct SortedNeighbours {
    std::vector<NeighbourPair> ids;
    std::vector<Vec3> dist;
};

// Keeps track of neighbour pair and
// corresponding kartesian distance
struct UnsortedNeighbour {

    // Pair of owner and corresponding neighbour ids
    NeighbourPair ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    Vec3 dist;
};

void neighbour_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const Scalar maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret);

void owner_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const Scalar maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret);

// A search cube stores first and last particle ids
struct SearchCube {
    size_t first;
    size_t last;
};

SortedNeighbours createNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Vec3> &pos,
    std::vector<SearchCube> &searchCubes);

struct SortedParticles {
    std::vector<SearchCube> searchCubes;
    std::vector<size_t> sorting_idxs;
    std::vector<Vec3> particles;
};

std::ostream &operator<<(std::ostream &os, NeighbourPair const &n);

std::vector<size_t> upper_neighbour_cubes(
    const SubDivision sub,
    const NeighbourIdHalfStencil &upper_neighbourId_stencil,
    const size_t id);

SortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Vec3> &unsorted_particles);

template <class A> class Field;

template <class A> class FieldAB;

using NeighbourFieldAB = FieldAB<Field<std::vector<NeighbourPair>>>;
#endif
