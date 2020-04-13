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

#include "Helper.hpp"
// #include <stdlib.h>     /* abort, NULL */

std::vector<Vec3> create_uniform_particle_plane(size_t n_particles) {
    std::vector<Vec3> points;
    points.reserve(n_particles);
    for (size_t j = 0; j < n_particles; j++) {
        for (size_t i = 0; i < n_particles; i++) {
            Scalar x = ((Scalar)i) / ((Scalar)n_particles);
            Scalar y = ((Scalar)j) / ((Scalar)n_particles);
            Scalar z = 0.0;
            points.push_back(Vec3 {x, y, z});
        }
    }
    return points;
}

std::vector<Vec3> create_uniform_particle_cube(
    Vec3 dimensions, Vec3 position, Scalar dx, Scalar noise) {

    srand(time(NULL));

    size_t nx {(size_t)((dimensions[0] - position[0]) / dx)};
    size_t ny {(size_t)((dimensions[1] - position[1]) / dx)};
    size_t nz {(size_t)((dimensions[2] - position[2]) / dx)};

    size_t ntot = nx * ny * nz;
    Scalar dx_layer = 0;

    std::vector<Vec3> points;
    points.reserve(ntot);
    for (size_t k = 0; k < nz; k++) {
        // on every second layer particles are moved dx/2 in x and y
        if (k % 2 == 0)
            dx_layer = 0.0;
        else
            dx_layer = dx / 2.0;

        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                Scalar nx = ((Scalar)(rand() % 100)) / 50.0 - 1.0;
                Scalar ny = ((Scalar)(rand() % 100)) / 50.0 - 1.0;
                Scalar nz = ((Scalar)(rand() % 100)) / 50.0 - 1.0;

                Scalar x = ((Scalar)i) * dx + nx * noise * dx + dx_layer;
                Scalar y = ((Scalar)j) * dx + ny * noise * dx + dx_layer;
                Scalar z = ((Scalar)k) * dx + nz * noise * dx;

                points.push_back(Vec3 {x, y, z});
            }
        }
    }
    return points;
}

std::vector<Vec3> create_uniform_particle_cube(size_t n_particles) {
    std::vector<Vec3> points;
    points.reserve(n_particles);
    for (size_t k = 0; k < n_particles; k++) {
        for (size_t j = 0; j < n_particles; j++) {
            for (size_t i = 0; i < n_particles; i++) {
                Scalar x = ((Scalar)i) / ((Scalar)n_particles);
                Scalar y = ((Scalar)j) / ((Scalar)n_particles);
                Scalar z = ((Scalar)k) / ((Scalar)n_particles);
                points.push_back(Vec3 {x, y, z});
            }
        }
    }
    return points;
}

Scalar rand01() { return ((Scalar)rand() / (RAND_MAX)); }

std::vector<Vec3> disperse_particles(std::vector<Vec3> &points, Scalar dx) {
    for (size_t piter = 0; piter < points.size(); piter++) {

        Scalar x = points[piter][0] + rand01() * dx;
        Scalar y = points[piter][1] + rand01() * dx;
        Scalar z = points[piter][2] + rand01() * dx;
        points[piter] = Vec3 {x, y, z};
    }
    return points;
}
