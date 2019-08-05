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

#ifndef HELPER_H
#define HELPER_H

#include <vector>     /* abort, NULL */
#include "CGALTYPEDEFS.hpp"


std::vector<Point> create_uniform_particle_plane(size_t n_particles);

std::vector<Point> create_uniform_particle_cube(size_t n_particles);

struct FixedDistanceParticles {

  std::vector<Point>  points; // id of original particle
  std::vector<size_t> fixId; // id of original particle
  std::vector<float>  maxDx; // max allowed distance

  // defines which motions are allowed
  // 0 - fixed, 1 - along line, 2 - plane, 3 half sphere,
  // 4 free
  std::vector<int>    mType;
  // motion vector for particle of mType
  // 0 - ignored, 1 direction of line, 2 normal to plane,
  // 3 - direction of half sphere, 4 0 -ignored

  std::vector<CGALVector>     dir;

  // Facet to which particle belongs to
  std::vector<Facet_handle>         facets;
};

float rand01();

std::vector<Point> disperse_particles(
        std::vector<Point>& points,
        float dx);

#endif
