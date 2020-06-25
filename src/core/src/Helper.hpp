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

#include <stddef.h> // for size_t
#include <vector>   // for vector

#include "Scalar.hpp"
#include "Vec3.hpp"

std::vector<Vec3> create_uniform_particle_plane(size_t n_particles);

std::vector<Vec3> create_uniform_particle_cube(size_t n_particles);

std::vector<Vec3> create_uniform_particle_cube(
    Vec3 dimensions, Vec3 position, Scalar dx, Scalar noise = 0);

Scalar rand01();

std::vector<Vec3> disperse_particles(std::vector<Vec3> &points, Scalar dx);

#endif
