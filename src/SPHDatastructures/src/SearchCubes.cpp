#ifndef SEARCHCUBES_H
#define SEARCHCUBES_H

#include "SPHDatastructures.hpp"

// inplace functions
void owner_cube_search(
    std::vector<Point> &pos,
    size_t first,
    size_t last,
    float maxDistanceSqr,
    SortedNeighbours &ret) {
    for (size_t oid = first; oid < last; oid++) {
        // local oid counter
        size_t oidctr = oid - first;

        const Point &opos = pos[oid];

        // starts from first+oidctr since oidctr are already
        // tested against this nid (distance pairs)
        // +1 to avoid testing oid==nid
        for (size_t nid = oid + 1; nid < last; nid++) {

            const Point &npos = pos[nid];

            const float distanceSqr = squared_distance(opos, npos);

            if(distanceSqr == 0) {continue;}

            if (distanceSqr < maxDistanceSqr) {
                ret.ownId.push_back(oid);
                ret.neighId.push_back(nid);
            };
        }
    }
};

std::array<bool, 27> vector_inner_owner_cube_search(
    Point opos,
    std::vector<Point> &pos,
    size_t first,
    size_t last,
    float maxDistanceSqr) {

    std::array<bool, 27> mask{};

    for (size_t i = 0; i < last - first; i++) {
        mask[i] = false;
    }

    for (size_t nid = first; nid < last; nid++) {

        // const Point npos = pos[nid];

        const float distanceSqr = squared_distance(opos, pos[nid]);

        mask[nid - first] = distanceSqr < maxDistanceSqr;
    }

    return mask;
}

void vector_owner_cube_search(
    std::vector<Point> &pos,
    size_t first,
    size_t last,
    float maxDistanceSqr,
    SortedNeighbours &ret) {
    for (size_t oid = first; oid < last; oid++) {
        // local oid counter

        const Point &opos = pos[oid];

        // starts from first+oidctr since oidctr are already
        // tested against this nid (distance pairs)
        // +1 to avoid testing oid==nid

        auto mask = vector_inner_owner_cube_search(
            opos, pos, oid + 1, last, maxDistanceSqr);

        size_t nid = first;
        for (size_t i = 0; i < last - oid + 1; i++) {
            nid++;
            if (mask[i]) {
                ret.ownId.push_back(oid);
                ret.neighId.push_back(nid);
            };
        }
    }
    // std::cout << "vector_owner_cube_search "
    //           << "ret.ownId.size()" << ret.ownId.size()
    //           << "ret.neighId.size()" << ret.neighId.size()
    //           << std::endl;
};

void neighbour_cube_search(
    std::vector<Point> &pos,
    size_t first,
    size_t last,
    size_t first_nc,
    size_t last_nc,
    float maxDistanceSqr,
    SortedNeighbours &ret) {

    // Step 3.1. set pivot particle
    for (size_t oid = first; oid < last; oid++) {
        const Point &opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for (size_t nid = first_nc; nid < last_nc; nid++) {
            const Point &npos = pos[nid];

            const float distanceSqr = squared_distance(opos, npos);

            if (distanceSqr < maxDistanceSqr) {
                ret.ownId.push_back(oid);
                ret.neighId.push_back(nid);
            };
        }
    }
};

size_t position_to_cube_id(SearchCubeDomain scd, const Point &p) {
    // TODO test if speed up n_cubes are copied to a const size_t nx ...

    const size_t nx = (size_t)scd.n.nx;
    const size_t ny = (size_t)scd.n.ny;
    const size_t nz = (size_t)scd.n.nz;

    // use idx instead of dx
    const float idx = scd.idx;
    const float mx = scd.min.x;
    const float my = scd.min.y;
    const float mz = scd.min.z;

    const size_t i {(size_t)((p.x() - mx)*idx)};
    const size_t j {(size_t)((p.y() - my)*idx)};
    const size_t k {(size_t)((p.z() - mz)*idx)};

    return i + nx * j + nx * ny * k;
}

SubDivision id_to_i_j_k(const size_t id, const SubDivision sub) {
    // TODO the nz-1 limiter seems unnecessary
    const size_t nx = sub.nx;
    const size_t ny = sub.ny;
    const size_t nz = sub.nz;
    const size_t nxny = nx * ny;
    const unsigned int k{(unsigned int)min(nz - 1, id / nxny)};
    const unsigned int j{(unsigned int)min(ny - 1, (id - k * nxny) / nx)};
    const unsigned int i{(unsigned int)min(nx - 1, id - j * nx - k * nxny)};
    // std::cout
    //           << [DEBUG ]
    //           << "id_to_i_j_k id " << id
    //           << " nx " << nx
    //           << " ny " << ny
    //           << " nz " << nz
    //           << " i " << i
    //           << " j " << j
    //           << " k " << k
    //           << std::endl;
    return {i, j, k};
}

// get non boundary neighbours from neighbour id, false - boundary, true - inner
// std::vector<bool> masked_neighbourId_stencil(
std::vector<bool>
lower_searchCubes_neighbour_mask(const size_t id, const SubDivision sub) {
    const SubDivision ijk = id_to_i_j_k(id, sub);

    std::vector<bool> ret(13, true);

    // std::cout << "lower_searchCubes_neighbour_mask id " << id
    //           << " i " << ijk.nx
    //           << " j " << ijk.ny
    //           << " k " << ijk.nz
    //           << std::endl;

    // mask left
    if (ijk.nx == 0)
        for (int i : {0, 3, 6, 9, 12}) ret[i] = false;
    // mask front
    if (ijk.ny == 0)
        for (int i : {0, 1, 2, 9, 10, 11}) ret[i] = false;
    // mask bottom
    if (ijk.nz == 0)
        for (int i : {0, 1, 2, 3, 4, 5, 6, 7, 8}) ret[i] = false;

    return ret;
}

std::vector<bool>
upper_searchCubes_neighbour_mask(const size_t id, const SubDivision sub) {
    const SubDivision ijk = id_to_i_j_k(id, sub);

    std::vector<bool> ret(13, true);

    // std::cout << "upper_searchCubes_neighbour_mask id " << id
    //           << " i " << ijk.nx
    //           << " j " << ijk.ny
    //           << " k " << ijk.nz
    //           << std::endl;

    // mask right
    // compares index ijk to number of subcubes, hence -1 is needed atm
    if (ijk.nx == sub.nx - 1) {
        for (int i : {0, 3, 6, 9, 12}) ret[i] = false;
    }
    // mask back
    if (ijk.ny == sub.ny - 1) {
        for (int i : {1, 2, 3, 10, 11, 12}) ret[i] = false;
    }
    // if mask top
    if (ijk.nz == sub.nz - 1) {
        for (int i : {4, 5, 6, 7, 8, 9, 10, 11, 12}) ret[i] = false;
    }

    return ret;
}

std::vector<size_t> lower_neighbour_cubes(
    /* const SearchCubeDomain& scd, */
    const SubDivision sub,
    const NeighbourIdHalfStencil &lower_neighbourId_stencil,
    const size_t id) {
    const std::vector<bool> mask = lower_searchCubes_neighbour_mask(id, sub);

    std::vector<size_t> neighbourIds;
    neighbourIds.reserve(13);

    // compute lower cube neighbour ids first
    for (size_t i = 0; i < lower_neighbourId_stencil.stencil.size(); i++) {

        // std::cout <<  "lower_neighbour_cubes"
        //           << " sub.nx " << sub.nx
        //           << " sub.ny " << sub.ny
        //           << " sub.nz " << sub.nz
        //           << " mask[i] " << mask[i]
        //           << " stencil[i] " <<
        //           lower_neighbourId_stencil.stencil[12-i]
        //           << std::endl;

        if (mask[i]) {
            const size_t nid {
                id - (size_t)lower_neighbourId_stencil.stencil[12 - i]};
            neighbourIds.push_back(nid);
        }
    }
    return neighbourIds;
}

std::vector<size_t> upper_neighbour_cubes(
    const SubDivision sub,
    const NeighbourIdHalfStencil &upper_neighbourId_stencil,
    const size_t id) {
    const std::vector<bool> mask = upper_searchCubes_neighbour_mask(id, sub);

    std::vector<size_t> neighbourIds;
    neighbourIds.reserve(13);

    for (size_t i = 0; i < upper_neighbourId_stencil.stencil.size(); i++) {
        if (mask[i]) {
            const size_t nid {(size_t)upper_neighbourId_stencil.stencil[i] +
                              id};
            neighbourIds.push_back(nid);
        }
    }
    return neighbourIds;
}

SortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Point> &unsorted_particles) {

  std::cout
    << "[DEBUG] SearchCubes.cpp 0 " << scd.nt
    << std::endl;
    std::vector<size_t> count(scd.nt, 0);

    const size_t n_particles = unsorted_particles.size();

    // get search cube id from particle position
    // the search cube id serves as integer for sorting

    // Step 1. setting counts
    for (size_t i = 0; i < n_particles; i++) {
        count[position_to_cube_id(scd, unsorted_particles[i])]++;
    }

    // Step 2. setting start index
    // this reuses the count array by summing up to a
    // continuous count
    std::vector<SearchCube> retc;
    retc.reserve(scd.nt);

    std::vector<size_t> rets;
    rets.reserve(n_particles);

    size_t first = 0;
    size_t last = 0;
    for (size_t i = 0; i < scd.nt; i++) {
        last = count[i];
        count[i] = first;
        retc.push_back({first, first + last});
        first += last;
    }

    std::vector<Point> retp(n_particles);
    for (size_t uid = 0; uid < n_particles; uid++) {
        // copy particle to new position in sorted array
        const size_t sid = position_to_cube_id(scd, unsorted_particles[uid]);

        // std::cout
        //   << "[DEBUG] SearchCubes.cpp::countingSortParticles"
        //   << " uid " << uid
        //   << " sid " << sid
        //   << " x " << unsorted_particles[uid].x()
        //   << " y " << unsorted_particles[uid].y()
        //   << " z " << unsorted_particles[uid].z()
        //   << std::endl;

        const size_t oid = count[sid]; // next

        rets.push_back(oid);

        retp[oid] = Point(unsorted_particles[uid]);
        count[sid]++;
    }

    return {retc, rets, retp};
};

SortedNeighbours mergedCountingSortAndNeighbourSearch(
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
  std::vector<SearchCube> retc;
  retc.reserve(scd.nt);

  size_t first = 0;
  size_t last = 0;
  for(size_t i=0; i<scd.nt; i++) {
    last = count[i];
    count[i] = first;

    retc.push_back({first, first+last});

    first += last;
  }

  SortedNeighbours ret {};
  ret.ownId.reserve(unsorted_particles.size());
  ret.neighId.reserve(40*unsorted_particles.size());

 // Iterate search cubes
  const float maxDistanceSqr = scd.dx*scd.dx;
  std::vector<Point> retp(n_particles);
  NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
  SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

  size_t last_sid = 0;
  for(size_t uid=0; uid<n_particles; uid++) {
    // copy particle to new position in sorted array
    const size_t sid = position_to_cube_id(scd, unsorted_particles[uid]);
    size_t oid = count[sid]; // next
    retp[oid] = Point(unsorted_particles[uid]);
    const Point& opos = retp[oid];

    const auto ncids = lower_neighbour_cubes(sub, ncidsten, sid);

    // search for neighbours in lower search cube neighbours
    for (size_t ncid: ncids) {

      const size_t first_nc = retc[ncid].first;
      // search only up to count
      const size_t last_nc = count[ncid];

      // Step 3.3. set iterate neighbour cube particle
      for (size_t nid=first_nc; nid<last_nc; nid++) {

        const Point &npos = retp[nid];

        const float distanceSqr =
          squared_distance(opos, npos);

        if (distanceSqr < maxDistanceSqr) {
          ret.ownId.push_back(oid);
          ret.neighId.push_back(nid);
        };
      }
    }

    // Wenn sid sid < last sid particle backfill neighbour cube forward search
    if (sid < last_sid) {
      // TODO stimmen die uncids bei backfill ??
      const auto uncids = upper_neighbour_cubes(sub, ncidsten, sid);

      // search for neighbours in lower search cube neighbours
      for (size_t uncid: ncids) {

        const size_t first_unc = retc[uncid].first;
        // search only up to count
        const size_t last_unc = count[uncid];

        // Step 3.3. set iterate neighbour cube particle
        for (size_t nid=first_unc; nid<last_unc; nid++) {

          const Point &npos = retp[nid];

          const float distanceSqr =
            squared_distance(opos, npos);

          if (distanceSqr < maxDistanceSqr) {
            ret.ownId.push_back(oid);
            ret.neighId.push_back(nid);
          };
        }
      }
    }

    last_sid = sid;

    // Wenn sid voll dann owner cube search
    if (count[sid]==retc[sid].last) {
      owner_cube_search(retp, first, last, maxDistanceSqr, ret);
    };

    count[sid]++;


    // if oid was last particle in current sid search all neighbour
    // pairs in sid

    // find pairs in current sid
    // iterate over particles in sid with lower oid
    // size_t first = retc[sid].first;
    // for (size_t nid=oid; nid>first; nid--) {
    //   const Point &npos = retp[nid];

    //   const float distanceSqr =
    //     squared_distance(opos, npos);

    //   if (distanceSqr < maxDistanceSqr) {
    //     ret.ownId.push_back(oid);
    //     ret.neighId.push_back(nid);
    //   };
    // }
  }

  return ret;
};

// void owner_cube_search(first, last, ret,) {};

SortedNeighbours
createNeighbours(SearchCubeDomain scd, SortedParticles &particles) {
    // Step 0 initialise return values

    SortedNeighbours ret {};
    /* ret.firstLast.reserve(particles.particles.size()); */
    ret.neighId.reserve(40 * particles.particles.size());
    ret.ownId.reserve(40 * particles.particles.size());

    std::vector<SearchCube> &searchCubes {particles.searchCubes};
    std::vector<Point> &pos {particles.particles};

    const float maxDistanceSqr = scd.dx * scd.dx;
    const size_t particle_stop_id = pos.size();

    NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
    SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

    // Step 1. get parent search cube
    for (size_t sid = 0; sid < scd.nt; sid++) {

        // Step 2. test all particle in parent search cube
        // Step 2.1 set first particle of parent cube
        const size_t first = searchCubes[sid].first;
        const size_t last = searchCubes[sid].last;

        owner_cube_search(pos, first, last, maxDistanceSqr, ret);

        const auto ncids = upper_neighbour_cubes(sub, ncidsten, sid);

        // Step 3. test all particles in neighbour search_cubes
        // Step 3.1. set neighbour search cube
        for (size_t ncid : ncids) {

            const size_t first_nc = searchCubes[ncid].first;
            const size_t last_nc = searchCubes[ncid].last;

            neighbour_cube_search(
                pos, first, last, first_nc, last_nc, maxDistanceSqr, ret);
        }
    }

    return ret;
};

    // Wrapper Class

    // class SearchCubes5 {

    // public:
    //   SearchCubes5(
    //                Logger logger,
    //                const std::vector<Point> &points,
    //                const float dx,
    //                bool keep_empty_cubes = true) {
    //     size_t n_points = std::cbrt(points.size());
    //     SearchCubeDomain scd = initSearchCubeDomain(n_points);
    //     SortedNeighbours sp = countingSortParticles(scd, points);
    //     // createNeighbours(scd, sp);
    //   };

    //   void sortParticles(const std::vector<Point> & unsorted_points) {
    //     size_t n_points = std::cbrt(unsorted_points.size());
    //     SearchCubeDomain scd = initSearchCubeDomain(n_points);
    //     SortedNeighbours sp = countingSortParticles(scd, unsorted_points);
    //   };
    // };

#endif
