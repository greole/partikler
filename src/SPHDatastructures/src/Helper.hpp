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
