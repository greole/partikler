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

// #include "Field.hpp" // for FieldAB
#include "cgal/CGALHelper.hpp"   // for STLSurfaceDist, Point3D
#include "cgal/CGALTYPEDEFS.hpp" // for Point, Facet_handle

#include <algorithm> // for std::sort

#include "Vec3.hpp"

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

struct SubDivision {
    // 3*4bytes = 12bytes
    // Max 65535**3 cubes, which are approx (65535*3)**3 particles
    unsigned int nx, ny, nz;
};

struct SearchCubeDomain {

    Point3D min;
    Point3D max;

    float dx;
    float idx;

    // 32 byte

    SubDivision n;
    size_t nt;

    // 24 byte

    /* size_t padding; */
};

SearchCubeDomain
initSearchCubeDomain(const std::vector<Point> &particles, float dx);

SearchCubeDomain
initSearchCubeDomain(const std::vector<Vec3> &particles, float dx);

struct NeighbourPair {
    size_t ownId;
    size_t neighId;
};

struct STLSortedNeighbours {
    std::vector<NeighbourPair> ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    std::vector<STLSurfaceDist> dist;
};

struct SortedNeighbours {
    std::vector<NeighbourPair> ids;
    std::vector<Vec3> dist;
};

struct STLUnsortedNeighbour {

    NeighbourPair ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    STLSurfaceDist dist;
};

struct UnsortedNeighbour {

    NeighbourPair ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    Vec3 dist;
};

void stl_neighbour_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const float maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret);

void stl_owner_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const float maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret);

void neighbour_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const float maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret);

void owner_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const float maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret);

// A search cube stores first and last particle ids
struct SearchCube {
    size_t first;
    size_t last;
};

STLSortedNeighbours createSTLNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> &pos,
    std::vector<SearchCube> &searchCubes,
    const std::vector<Facet_handle> &facets);

SortedNeighbours createNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Vec3> &pos,
    std::vector<SearchCube> &searchCubes);

// SortedNeighbours createNeighbours(
//     const SearchCubeDomain scd,
//     const std::vector<Point> &pos,
//     std::vector<SearchCube> &searchCubes);

struct STLSortedParticles {
    std::vector<SearchCube> searchCubes;
    std::vector<size_t> sorting_idxs;
    std::vector<Point> particles;
};

struct SortedParticles {
    std::vector<SearchCube> searchCubes;
    std::vector<size_t> sorting_idxs;
    std::vector<Vec3> particles;
};

STLSortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Point> &unsorted_particles);

SortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Vec3> &unsorted_particles);

template <class A> class Field;

template <class A> class FieldAB;

using NeighbourFieldAB = FieldAB<Field<std::vector<NeighbourPair>>>;
#endif
