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

#ifndef PARTIKLER_CGALSEARCHCUBES_INCLUDED
#define PARTIKLER_CGALSEARCHCUBES_INCLUDED

#include <omp.h>    // omp_get_num_threads
#include <stddef.h> // for size_t
#include <vector>   // for vector

// #include "Field.hpp" // for FieldAB
#include "cgal/CGALHelper.hpp"   // for STLSurfaceDist, Point3D
#include "cgal/CGALTYPEDEFS.hpp" // for Point, Facet_handle

#include <algorithm> // for std::sort

#include "Scalar.hpp"
#include "SearchCubes.hpp"
#include "Vec3.hpp"

struct STLSortedNeighbours {
    std::vector<NeighbourPair> ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    std::vector<STLSurfaceDist> dist;
};

struct STLUnsortedNeighbour {

    NeighbourPair ids;

    // since computation of neighbouring particles
    // on different STL is expensive a vector holding
    // the STLSurfaceDist is stored here
    STLSurfaceDist dist;
};

std::pair<Vec3, Vec3> cgal_bounding_box(std::vector<Point> const &particles);

void stl_neighbour_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const Scalar maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret);

void stl_owner_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const Scalar maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret);

STLSortedNeighbours createSTLNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> &pos,
    std::vector<SearchCube> &searchCubes,
    const std::vector<Facet_handle> &facets);

struct STLSortedParticles {
    std::vector<SearchCube> searchCubes;
    std::vector<size_t> sorting_idxs;
    std::vector<Point> particles;
};

STLSortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Point> &unsorted_particles);

template <class A> class Field;

template <class A> class FieldAB;

using NeighbourFieldAB = FieldAB<Field<std::vector<NeighbourPair>>>;
#endif
