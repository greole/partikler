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

#include "SPHDatastructures.hpp"
#include <algorithm>    // std::min
#include <math.h>    // std::min
#include <time.h>    // std::min
#include <stdlib.h>    // std::min

void compute_forces(
    Logger logger,
    const searchcubes::SortedNeighbours &particle_neighbours,
    const SPHPointField &particles,
    const std::vector<Facet_handle> & facets,
    const float dx,
    SPHVectorField &f)
{
  // TODO runtime.submodels("force").calculate()
  // Lennard Jones Potential
  logger.info_begin() << "Computing forces";
  f.set(zeroVec);

  const size_t size {particle_neighbours.ids.size()};
  const float sigma = dx;
  const float epsilon = 1e-06;

  for (size_t pid = 0; pid < size; pid++) {

    size_t oid = particle_neighbours.ids[pid].ownId;
    size_t nid = particle_neighbours.ids[pid].neighId;
    const Point &opos = particles[oid];
    const Point &npos = particles[nid];

    Facet F1 = Facet(*facets[oid]);
    Facet F2 = Facet(*facets[nid]);

    CGALVector lenVo = particle_neighbours.dist[pid].on;
    CGALVector lenVn = particle_neighbours.dist[pid].no;

    double len = (double) particle_neighbours.dist[pid].len;

    double rat = min(4.0, max(0.01, sigma/len));
    double fi = 4.*epsilon*(pow(rat, 12.) - pow(rat, 6.) );

      for(int j=0;j<3;j++) {
        if (lenVo[j] != 0.0) {
          f[oid][j] -= (float) fi/lenVo[j];
        }
        if (lenVn[j] != 0.0) {
          f[nid][j] -= (float) fi/lenVn[j];
        }
    }
  }
  logger.info_end();
};

void add_random_noise(Logger logger, SPHVectorField &u, const float dx) {
}

void limit_dt_du(
         const SPHVectorField &du,
         const float maxDx,
         float & dt
         ) {
  const float maxCFL = 0.5; // DONT HARDCODE
  float maxDu = du.norm().get_max();
  std::cout << "maxDu " << maxDu << std::endl;
  float CFL = maxDu*dt*dt/maxDx;
  std::cout << "maxCFL " << CFL << std::endl;
  // max double timestep
  float two = 2.0;
  float change = min(two, maxCFL/CFL);
  dt = dt * change;
}

float limit_sphere_frac(CGALVector OP, CGALVector PN, float maxDx2) {
    float a = PN.squared_length();
    float b = 2.0 * OP * PN;
    float c = OP.squared_length() - maxDx2;

    float lambda1 = (-b + std::sqrt(b * b - 4 * a * c)) / (2.0 * a);
    float unity = 1.0;
    float zero = 0.0;
    return std::max(zero, std::min(unity, lambda1));
}

CGALVector
half_sphere(Point O, Point P, CGALVector PN, CGALVector D, float maxDx) {
    Point N = P + PN;
    CGALVector ON = N - O;

    float maxDx2 = maxDx * maxDx;

    // if new particle position is inside half sphere everything is fine
    // Positive dot product means new particle position is
    // in upper half of the sphere, as defined by extrusion dir
    bool newPosIsUpperHalfSphere = (ON * D >= 0);
    // Distance to original point is less than maxDx
    bool newPosIsInsideSphere = ON.squared_length() <= maxDx2;
    if (newPosIsUpperHalfSphere && newPosIsInsideSphere) {
        return PN;
    };

    CGALVector OP = P - O;
    float a = PN.squared_length();
    bool currentPosInsideSphere = OP.squared_length() <= maxDx2;

    // if particle is currently inside sphere but new position
    // would be outside the sphere dx is to be limited
    float frac = 0;
    if (currentPosInsideSphere && (a > 0)) {
        frac = limit_sphere_frac(OP, PN, maxDx2);

        // if dot product of extrusion dir and dx is negative
        // particle travels towards half sphere boundary thus frac
        // is potentially to be limited
        float denom = PN * D;
        bool canCrossHalfSphereBound = denom < 0;
        if (canCrossHalfSphereBound) {
            float nom = -OP * D;
            float lambda2 = nom / denom;
            float zero = 0.0;
            frac = std::max(zero, std::min(frac, lambda2));
        }
    }
    return PN * frac;
}


void snap_to_surface(
  std::vector<Facet_handle>& facets,
  SPHPointField& pos
                     ){

  std::vector<Vector> corrector(pos.size(), {0,0,0});

  for (size_t i = 0; i < pos.size(); i++) {
  // compute distance to facet
    // Facet facet = Facet(*facets[i]);
    // CGALVector N = facet_normal(facet);
    // float d = facet.plane().d();
    // CGALVector P = pos[i]-Point(0, 0, 0);
    // float dist = (N*P+d)/std::sqrt(N.squared_length());

    Facet facet = Facet(*facets[i]);
    auto h = facet.facet_begin();

    Plane p {
             h->vertex()->point(),
             h->next()->vertex()->point(),
             h->next()->next()->vertex()->point(),
    };

    Point P = p.projection(pos[i]);
    CGALVector dx = P - pos[i];

    // if (dist != 0.0) {
    //   std::cout
    //     << "[DEBUG SLIDE] Snap to STL"
    //     << " Particle Id " << i
    //     << " Pos " << P
    //     << " dist " << dist
    //     << " N " << N
    //     << " d*N" << d*N
    //     << " p-d*N " << pos[i] - d*N
    //     << " DONE "
    //     << std::endl;
        corrector[i] = Vector{(float)dx.x(), (float)dx.y(), (float)dx.z()};
  }

  SPHVectorField cor = SPHVectorField(corrector, {"tmpx", "tmpy", "tmpz"}, "tmp");
}
