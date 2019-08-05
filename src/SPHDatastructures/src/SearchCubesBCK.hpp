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

#ifndef SEARCHCUBE5_H
#define SEARCHCUBE5_H

struct SearchCubeDomain {

  float x_min;
  float y_min;
  float z_min;

  float x_max;
  float y_max;
  float z_max;

  float dx;
  float idx;

  // 32 byte

  size_t nx;
  size_t ny;
  size_t nz;
  size_t nt;

  // 24 byte

  /* size_t padding; */

};

SearchCubeDomain initSearchCubeDomain(const size_t n_particles) {

  const size_t ncubes = std::max((size_t)1, (size_t)ceil((float)n_particles/3.0));
  const float idx = (float)ncubes;
  const float dx = 1.0 /idx;
  return SearchCubeDomain {
    0.0,
    0.0,
    0.0,
    1.0,
    1.0,
    1.0,
    dx,
    idx,
    ncubes, ncubes, ncubes,
    ncubes*ncubes*ncubes,
      };
}

size_t position_to_cube_id(SearchCubeDomain scd, const Point& p) {
  // TODO scale bounding box a bit
  // TODO test if speed up n_cubes are copied to a const size_t nx ...

  const size_t nx = scd.nx;
  const size_t ny = scd.ny;
  const size_t nz = scd.nz;

  // use idx instead of dx
  const float idx = scd.idx;
  const float mx = scd.x_min;
  const float my = scd.y_min;
  const float mz = scd.z_min;

  const size_t i {std::min(nx, (size_t)((p.x() - mx)*idx))};
  const size_t j {std::min(ny, (size_t)((p.y() - my)*idx))};
  const size_t k {std::min(nz, (size_t)((p.z() - mz)*idx))};

  /* const size_t i {(size_t)((p.x() - mx)*idx)}; */
  /* const size_t j {(size_t)((p.y() - my)*idx)}; */
  /* const size_t k {(size_t)((p.z() - mz)*idx)}; */

  return i + nx*j + nx*ny*k;
}


std::vector<size_t> id_to_i_j_k (
   const size_t id, size_t nx, size_t ny, size_t nz)
{
  // TODO the nz-1 limiter seems unnecessary
  const size_t k {min( nz - 1, id / nx*ny)};
  const size_t j {min( ny - 1, id / nx)};
  const size_t i {min( nx - 1, id - j*nx - k*nx*ny)};
  return {i, j, k};
}

// get non boundary neighbours from neighbour id, false - boundary, true - inner
std::vector<bool> masked_neighbourId_stencil(
  const size_t id, size_t nx, size_t ny, size_t nz
) {
    const std::vector<size_t> ijk = id_to_i_j_k(
      id, nx, ny, nz
    );
    std::vector<bool> ret (13, true);

    if  (ijk[0] == nx - 1) {
      for (int i : {0, 3, 6, 9, 12}) {
        ret[i] = false;
      }
    }
    // if j == 0 mask back
    if  (ijk[1] == ny - 1) {
      for (int i : {1, 2, 3, 10, 11, 12}) {
        ret[i] = false;
      }
    }
    // if k == 0 mask top
    if  (ijk[2] == nz - 1) {
      for (int i : {4, 5, 6, 7, 8, 9, 10, 11}) {
        ret[i] = false;
      }
    }
    return ret;
}

struct NeighbourIdStencil {
  // computes the symmetric upper part of the neighbour id stencil

  // In generall neighbour ids are as following
  //  n_cubes_[1]=3
  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |06 07 08     |15 16 17     |24 25 26 // back
  //  |03 04 05     |12 13 14     |21 22 23 // middle
  //  |00 01 02     |09 10 11     |18 19 20 // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3
  //
  //  Bottom        Middle        Top
  //  ^ y           ^ y           ^ y
  //  |06 07 08     |1  2  3      |10 11  12 // back
  //  |03 04 05     |x  x  0      |7   8   9 // middle
  //  |00 01 02     |x  x  x      |4   5   6 // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3
  // search cube ids start from min(x) to max(z)
  // in a rhs coordinate system
  // - the right neighbour (x-axis) is +1
  // - the back neighbour  (y-axis) is +(n_cubes_[0])
  // - the upper neighbour (z-axis) is +(n_cubes_[0]*n_cubes[1])
  std::vector<size_t> stencil;

  NeighbourIdStencil (size_t nx, size_t ny, size_t nz) {
    // TODO leave it a size_t, iterate only till 12, since
    // the stencil is symmetric anyway
    const size_t ny_nx = nx*ny;

    stencil = std::vector<size_t> {
                      1, // right
                 nx - 1, // back
                 nx    , // back
                 nx + 1, // back
         ny_nx - nx - 1, // back
         ny_nx - nx    , // back
         ny_nx - nx + 1, // back
         ny_nx       -1, // back
         ny_nx         , // back
         ny_nx      + 1, // back
         ny_nx + nx - 1, // back
         ny_nx + nx    , // back
         ny_nx + nx + 1  // back
    };
  };
};

std::vector<size_t> neighbour_cubes (
      /* const SearchCubeDomain& scd, */
      const size_t nx,
      const size_t ny,
      const size_t nz,
      const NeighbourIdStencil& neighbourId_stencil,
      const size_t id
) {
  const std::vector<bool> neighbour_mask = masked_neighbourId_stencil(
    id, nx, ny, nz
    );

    std::vector<size_t> neighbourIds;
    neighbourIds.reserve(14);

    // compute lower cube neighbour ids first
    /* for (size_t i = 0; i<neighbourId_stencil.stencil.size(); i++) */
    /* { */
    /*     if (neighbour_mask[12-i]) { */
    /*         // TODO remove this check once neighbour_mask works properly */
    /*         // check if nid underflows size_t */
    /*         if (neighbourId_stencil.stencil[i] <= id) { */
    /*             const size_t nid {id - neighbourId_stencil.stencil[i]}; */
    /*             neighbourIds.push_back(nid); */
    /*         } */
    /*     } */
    /* } */

    // compute upper cube neighbour ids second
    for (size_t i = 0; i<neighbourId_stencil.stencil.size(); i++)
    {
        if (neighbour_mask[i]) {
            const size_t nid {neighbourId_stencil.stencil[i] + id};
            neighbourIds.push_back(nid);
        }
    }
    return neighbourIds;
}


struct SearchCube5 {

  size_t first;
  size_t last;

  SearchCube5(size_t first, size_t last) : first(first), last(last){};
};

struct SortedParticles {
  std::vector<SearchCube5> searchCubes;
  std::vector<Point> particles;
};

SortedParticles countingSortParticles(
      const SearchCubeDomain scd,
      const std::vector<Point>& unsorted_particles) {

  std::vector<int> count (scd.nt, 0);
  const size_t n_particles = unsorted_particles.size();

  // get search cube id from particle position
  // the search cube id serves as integer for sorting

  // Step 1. setting counts
  for(size_t i=0; i<n_particles; i++) {
    count[position_to_cube_id(scd, unsorted_particles[i])]++;
  }

  // Step 2. setting start index
  // this reuses the count array by summing up to a
  // continuous count
  std::vector<SearchCube5> retc;
  retc.reserve(scd.nt);

  size_t first = 0;
  size_t last = 0;
  for(size_t i=0; i<scd.nt; i++) {
    last = count[i];
    count[i] = first;

    retc.push_back({first, first+last});

    first += last;
  }

  // Step 3. Writing output array
  std::vector<Point> retp(n_particles);
  for(size_t i=0; i<n_particles; i++) {
    size_t next = count[position_to_cube_id(scd, unsorted_particles[i])];
    retp[next] = Point(unsorted_particles[i]);
    count[position_to_cube_id(scd, unsorted_particles[i])]++;
  }

  return SortedParticles {retc, retp};
};

struct SortedNeighbours {
  // Each particle id has a first and last neighbour
  // stored in firstLast, here SearchCube5 is reused
  /* std::vector<packed> ownId; */
  std::vector<size_t> ownId;
  std::vector<size_t> neighId;
};

struct SortedNeighboursNElems {
  // Each particle id has a first and last neighbour
  // stored in firstLast, here SearchCube5 is reused
  std::vector<int> nElems;
  std::vector<size_t> neighbId;
};

SortedNeighbours createNeighbours(
      SearchCubeDomain scd,
      SortedParticles &particles
                                    ) {
  // Step 0 initialise return values

  SortedNeighbours ret {};
  /* ret.firstLast.reserve(particles.particles.size()); */
  ret.neighId.reserve(40*particles.particles.size());
  ret.ownId.reserve(40*particles.particles.size());

  std::vector<SearchCube5>& searchCubes {particles.searchCubes};
  std::vector<Point>& pos {particles.particles};

  const float maxDistanceSqr = scd.dx*scd.dx;
  const size_t particle_stop_id = pos.size();

  NeighbourIdStencil ncidsten(scd.nx, scd.ny, scd.nz);

  // Step 1. get parent search cube
  for (size_t sid=0; sid<scd.nt; sid++) {

    const auto ncids = neighbour_cubes(scd.nx, scd.ny, scd.nz, ncidsten, sid);

    // Step 2. test all particle in parent search cube
    // Step 2.1 set first particle of parent cube 
    const size_t first = searchCubes[sid].first;
    const size_t last = searchCubes[sid].last;
    // number of particles in search cube

    for  (size_t oid=first; oid<last; oid++) {

      // local oid counter
      size_t oidctr = oid-first;

      const Point opos = pos[oid];

      // starts from first+oidctr since oidctr are already
      // tested against this nid (distance pairs)
      // +1 to avoid testing oid==nid
      for (size_t nid=first+oidctr+1; nid<last; nid++) {

        size_t nidctr = nid - first;

        const Point npos = pos[nid];

        const float distanceSqr =
            squared_distance(opos, npos);

        if (distanceSqr < maxDistanceSqr) {
          ret.ownId.push_back(oid);
          ret.neighId.push_back(nid);
        };
      }
    }

    // Step 3. test all particles in neighbour search_cubes
    // Step 3.1. set neighbour search cube
    for (size_t ncid: ncids) {

      const size_t first_nc = searchCubes[ncid].first;
      const size_t last_nc = searchCubes[ncid].last;

      // Step 3.1. set pivot particle
      for  (size_t oid=first; oid<last; oid++) {
        // size_t oidctr = oid-first;

        const Point opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for  (size_t nid=first_nc; nid<last_nc; nid++) {

          size_t nidctr = nid - first_nc;

          const Point npos = pos[nid];

          const float distanceSqr =
            squared_distance(opos, npos);

          if (distanceSqr < maxDistanceSqr) {
            ret.ownId.push_back(oid);
            ret.neighId.push_back(nid);
          };
        }
      }
    }
  }

  return ret;
};

SortedNeighboursNElems createNeighboursPivotParticle(
      SearchCubeDomain scd,
      SortedParticles &particles
                                    ) {
  // Step 0 initialise return values

  SortedNeighboursNElems ret {};
  ret.nElems.reserve(particles.particles.size());
  ret.neighbId.reserve(40*particles.particles.size());

  // Step 1. get parent search cube
  std::vector<SearchCube5>& searchCubes {particles.searchCubes};
  std::vector<Point>& pos {particles.particles};
  const float maxDistanceSqr = scd.dx*scd.dx;

  const size_t particle_stop_id = pos.size();

  NeighbourIdStencil ncidsten(scd.nx, scd.ny, scd.nz);

  size_t sid = 0;

  auto ncids = neighbour_cubes(scd.nx, scd.ny, scd.nz, ncidsten, sid);
  for (size_t oid=0; oid<particles.particles.size(); oid++) {

    // Step 2. test all particle in parent search cube
    // Step 2.1 set first particle of parent cube 
    const Point &opos = pos[oid];

    // check current parent search cube
    while (oid > searchCubes[sid].last) {
      sid++;
      // update ncids only if search cube changed
      ncids = neighbour_cubes(scd.nx, scd.ny, scd.nz, ncidsten, sid);
    }

    const size_t first = searchCubes[sid].first;
    const size_t last = searchCubes[sid].last;

    // starts from first+oidctr since oidctr are already
    // tested against this nid (distance pairs)
    // +1 to avoid testing oid==nid
    for (size_t nid=first; nid<last; nid++) {

      if (nid == oid) continue;

      const Point npos = pos[nid];

      const float distanceSqr =
          squared_distance(opos, npos);

      if (distanceSqr < maxDistanceSqr) {
        // store in temporary array until all neighbours for particle oid are
        // found, store nid and oid at once
        ret.neighbId.push_back(nid);
      };
    }

    // Step 3. test all particles in neighbour search_cubes
    // Step 3.1. set neighbour search cube
    /* for (int ncid: neighbour_cubes(scd, nids, sid) ) { */
    for (size_t ncid: ncids) {

      const size_t first_nc = searchCubes[ncid].first;
      const size_t last_nc = searchCubes[ncid].last;

      // Step 3.1. set pivot particle
      // Step 3.3. set iterate neighbour cube particle
      for  (size_t nid=first_nc; nid<last_nc; nid++) {

        const Point npos = pos[nid];

        const float distanceSqr =
          squared_distance(opos, npos);

        if (distanceSqr < maxDistanceSqr) {
          ret.neighbId.push_back(nid);
        };
      }
    }
  }
  return ret;
};

// Wrapper Class

class SearchCubes5 {

public:
  SearchCubes5(
               Logger logger,
               const std::vector<Point> &points,
               const float dx,
               bool keep_empty_cubes = true) {
    size_t n_points = std::cbrt(points.size());
    SearchCubeDomain scd = initSearchCubeDomain(n_points);
    SortedParticles sp = countingSortParticles(scd, points);
    createNeighbours(scd, sp);
  };

  void sortParticles(const std::vector<Point> & unsorted_points) {
    size_t n_points = std::cbrt(unsorted_points.size());
    SearchCubeDomain scd = initSearchCubeDomain(n_points);
    SortedParticles sp = countingSortParticles(scd, unsorted_points);
  };
};

#endif
