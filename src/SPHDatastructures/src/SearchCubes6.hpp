#ifndef SEARCHCUBE6_H
#define SEARCHCUBE6_H

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

SearchCubeDomain initSearchCubeDomain2(const size_t n_particles) {

  const size_t ncubes = std::max((size_t)1, (size_t)ceil((float)n_particles/3.0));
  const float idx = (float)ncubes;
  const float dx = 1.0/idx;
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
    // TODO test if wrapping it in a loop causes performance degradation
    // TODO can probably generated with constexpr

    // if i == 0 mask right
    if  (ijk[0] == nx - 1) {
        for (int i : {0, 1, 4, 7, 10, 13}) {
              ret[i] = false;
        }
    }
    // if j == 0 mask back
    if  (ijk[1] == ny - 1) {
        for (int i : {2, 3, 4, 11, 12, 13}) {
            ret[i] = false;
        }
    }
    // if k == 0 mask top
    if  (ijk[2] == nz - 1) {
        for (int i : {5, 6, 7, 8, 9, 10, 11, 12, 13}) {
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
  //  |06 07 08     |2  3  4      |11 12  13 // back
  //  |03 04 05     |x  x  1      |8   9  10 // middle
  //  |00 01 02     |x  x  0      |5   6   7 // front
  //  |-------> x   |-------> x   |-------> x
  //  left  right                           n_cubes_[0]=3
  //
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
    neighbourIds.reserve(13);

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


struct SearchCube6 {

  size_t first;
  size_t last;

  SearchCube6(size_t first, size_t last) : first(first), last(last){};
};

struct SortedParticles2 {
  std::vector<SearchCube6> searchCubes;
  std::vector<Point> particles;
};

SortedParticles2 countingSortParticles2(
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
  std::vector<SearchCube6> retc;
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

struct SortedNeighbours2 {
  // Each particle id has a first and last neighbour
  // stored in firstLast, here SearchCube6 is reused
  /* std::vector<SearchCube6> firstLast; */
  /* std::vector<size_t> neighbId; */
  std::vector<size_t> ownId;
  std::vector<size_t> neighId;
};

/* struct SortedNeighbours2NElems { */
/*   // Each particle id has a first and last neighbour */
/*   // stored in firstLast, here SearchCube6 is reused */
/*   /\* std::vector<int> nElems; *\/ */
/*   std::vector<size_t> ownerId; */
/*   std::vector<size_t> neighId; */
/* }; */

SortedNeighbours2 createNeighbours2(
      SearchCubeDomain scd,
      SortedParticles &particles
                                    ) {
  // Step 0 initialise return values

  SortedNeighbours2 ret {};
  ret.neighId.reserve(40*particles.particles.size());
  ret.ownId.reserve(40*particles.particles.size());

  std::vector<SearchCube6>& searchCubes {particles.searchCubes};
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

        const Point npos = pos[nid];

        const float distanceSqr =
            squared_distance(opos, npos);

        if (distanceSqr < maxDistanceSqr) {
          // store in temporary array until all neighbours for particle oid are
          // found, store nid and oid at once
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
        size_t oidctr = oid-first;

        const Point opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for  (size_t nid=first_nc+oidctr; nid<last_nc; nid++) {

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



class SearchCubes6 {

public:
  SearchCubes6(
               Logger logger,
               const std::vector<Point> &points,
               const float dx,
               bool keep_empty_cubes = true) {
    size_t n_points = std::cbrt(points.size());
    SearchCubeDomain scd = initSearchCubeDomain2(n_points);
    SortedParticles2 sp = countingSortParticles2(scd, points);
    createNeighbours2(scd, sp);
  };

  void sortParticles(const std::vector<Point> & unsorted_points) {
    size_t n_points = std::cbrt(unsorted_points.size());
    SearchCubeDomain scd = initSearchCubeDomain2(n_points);
    SortedParticles2 sp = countingSortParticles2(scd, unsorted_points);
  };
};

#endif
