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

std::vector<Point> create_uniform_particle_plane(size_t n_particles) {
    std::vector<Point> points;
    points.reserve(n_particles);
        for (size_t j = 0; j < n_particles; j++) {
            for (size_t i = 0; i < n_particles; i++) {
                float x = ((float)i) / ((float)n_particles);
                float y = ((float)j) / ((float)n_particles);
                float z = 0.0;
                points.push_back(Point(x, y, z));
            }
        }
    return points;
}

std::vector<Point> create_uniform_particle_cube(size_t n_particles) {
    std::vector<Point> points;
    points.reserve(n_particles);
    for (size_t k = 0; k < n_particles; k++) {
        for (size_t j = 0; j < n_particles; j++) {
            for (size_t i = 0; i < n_particles; i++) {
                float x = ((float)i) / ((float)n_particles);
                float y = ((float)j) / ((float)n_particles);
                float z = ((float)k) / ((float)n_particles);
                points.push_back(Point(x, y, z));
            }
        }
    }
    return points;
}

float rand01() { return ((float)rand() / (RAND_MAX)); }

std::vector<Point> disperse_particles(
        std::vector<Point>& points,
        float dx)
{
    for (size_t piter=0; piter<points.size(); piter++) {

        float x = points[piter].x() + rand01()*dx;
        float y = points[piter].x() + rand01()*dx;
        float z = points[piter].x() + rand01()*dx;
        points[piter] = Point(x,y,z);
    }
    return points;
}
