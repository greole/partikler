#ifndef HELPER_H
#define HELPER_H

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

#endif
