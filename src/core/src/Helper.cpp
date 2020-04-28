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

    Scalar xpos = 0;
    Scalar xlen = dimensions[0] - position[0];
    Scalar ypos = 0;
    Scalar ylen = dimensions[1] - position[1];
    Scalar zpos = 0;
    Scalar zlen = dimensions[2] - position[2];

    size_t nx {(size_t)((dimensions[0] - position[0]) / dx)};
    size_t ny {(size_t)((dimensions[1] - position[1]) / dx)};
    size_t nz {(size_t)((dimensions[2] - position[2]) / dx)};

    size_t ntot = nx * ny * nz;
    Scalar dx_layer = 0;

    std::vector<Vec3> points;
    points.reserve(ntot);

    // TODO this currently has fringes on its edges
    size_t j = 0;
    size_t k = 0;
    while (zpos < zlen) {

        if (j % 2 == 0)
            ypos = position[1];
        else
            ypos = position[1] + dx / 2.0;

        while (ypos < ylen) {

            if (k % 2 == 0)
                xpos = position[0];
            else
                xpos = position[0] + dx / 2.0;

            while (xpos < xlen) {
                points.push_back(Vec3 {xpos, ypos, zpos});
                xpos += dx;
            }

            ypos += dx/2.0;
            k++;
        }
        k = 0;
        zpos += dx/2.0;
        j++;
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
