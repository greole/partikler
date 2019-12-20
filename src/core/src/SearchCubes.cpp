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


#include "SearchCubes.hpp"


SearchCubeDomain initSearchCubeDomain(
    const std::vector<Point> & particles, float dx) {
    size_t n_particles = particles.size();
    auto bound_box = bounding_box(particles.begin(), particles.end());

    // bounds are scaled a bit to avoid particles on exact boundary
    const float domain_extrusion = 0.01 * dx; // in dx
    float min_x = bound_box.min().x() - domain_extrusion;
    float min_y = bound_box.min().y() - domain_extrusion;
    float min_z = bound_box.min().z() - domain_extrusion;
    float max_x = bound_box.max().x() + domain_extrusion;
    float max_y = bound_box.max().y() + domain_extrusion;
    float max_z = bound_box.max().z() + domain_extrusion;
    // std::cout << "initSearchCubeDomain" << " min_z " << min_z << std::endl;
    // std::cout << "initSearchCubeDomain" << " max_z " << max_z << std::endl;
    // for (auto p: particles) std::cout << p << std::endl;

    // ceil with fixed dx increases max(x,y,z) a bit
    unsigned int n_cubes_x = (unsigned int)ceil(((max_x - min_x) / dx));
    unsigned int n_cubes_y = (unsigned int)ceil(((max_y - min_y) / dx));
    unsigned int n_cubes_z = (unsigned int)ceil(((max_z - min_z) / dx));

    // Check if domain is not degenerated
    // search cubes currently assume 26 neighbour cubes
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_x >= 3)
        << " n_cubes_x less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_y >= 3)
        << " n_cubes_y less than 3";
    Logger(3, "initSearchCubeDomain:255").check(n_cubes_z >= 3)
        << " n_cubes_z less than 3";

    std::cout << " [DEBUG] SPHDatastructures.hpp::initSearchCubeDomain"
              << " min_x " << min_x << " min_y " << min_y << " min_z " << min_z
              << " max_x " << max_x << " max_y " << max_y << " max_z " << max_z
              << " n_cubes_x " << n_cubes_x << " n_cubes_y " << n_cubes_y
              << " n_cubes_z " << n_cubes_z << " dx " << dx << std::endl;

    const float idx = 1.0 / dx;
    return SearchCubeDomain {
        Point3D {min_x, min_y, min_z},
        Point3D {max_x, max_y, max_z},
        dx,
        idx,
        SubDivision {n_cubes_x, n_cubes_y, n_cubes_z},
        (size_t)n_cubes_x * (size_t)n_cubes_y * (size_t)n_cubes_z,
    };
}

// inplace functions
void owner_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const float maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    SortedNeighbours &ret) {

    for (size_t oid = first; oid < last; oid++) {
        const Point &opos = pos[oid];

        // starts from first+oidctr since oidctr are already
        // tested against this nid (distance pairs)
        // +1 to avoid testing oid==nid
        for (size_t nid = oid + 1; nid < last; nid++) {
            const Point &npos = pos[nid];

            if (opos == npos) {
                continue;
            }

            STLSurfaceDist sd =
                compute_STLSurface_dist(opos, npos, facets[oid], facets[nid]);

            const float distanceSqr = sd.len * sd.len;

            if (distanceSqr < maxDistanceSqr) {
                ret.ids.push_back({oid, nid});
                ret.dist.push_back(sd);
            };
        }
    }
}

std::array<bool, 27> vector_inner_owner_cube_search(
    Point opos,
    std::vector<Point> &pos,
    size_t first,
    size_t last,
    float maxDistanceSqr) {

    std::array<bool, 27> mask {};

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

void neighbour_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const float maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    SortedNeighbours &ret) {

    // Step 3.1. set pivot particle
    for (size_t oid = first; oid < last; oid++) {
        const Point &opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for (size_t nid = first_nc; nid < last_nc; nid++) {
            const Point &npos = pos[nid];

            STLSurfaceDist sd =
                compute_STLSurface_dist(opos, npos, facets[oid], facets[nid]);

            const float distanceSqr = sd.len * sd.len;

            if (distanceSqr < maxDistanceSqr) {
                ret.ids.push_back({oid, nid});
                ret.dist.push_back(sd);
            };
        }
    }
}

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

    const size_t i {(size_t)((p.x() - mx) * idx)};
    const size_t j {(size_t)((p.y() - my) * idx)};
    const size_t k {(size_t)((p.z() - mz) * idx)};

    return i + nx * j + nx * ny * k;
}

SubDivision id_to_i_j_k(const size_t id, const SubDivision sub) {
    // TODO the nz-1 limiter seems unnecessary
    const size_t nx = sub.nx;
    const size_t ny = sub.ny;
    const size_t nz = sub.nz;
    const size_t nxny = nx * ny;
    const unsigned int k {(unsigned int)min(nz - 1, id / nxny)};
    const unsigned int j {(unsigned int)min(ny - 1, (id - k * nxny) / nx)};
    const unsigned int i {(unsigned int)min(nx - 1, id - j * nx - k * nxny)};
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
    //           << " sub.nx " << sub.nx
    //           << " sub.ny " << sub.ny
    //           << " sub.nz " << sub.nz
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
}

SortedNeighbours mergedCountingSortAndNeighbourSearch(
    const SearchCubeDomain scd, const std::vector<Point> &unsorted_particles) {

    //  std::vector<int> count (scd.nt, 0);
    //  const size_t n_particles = unsorted_particles.size();

    //  // get search cube id from particle position
    //  // the search cube id serves as integer for sorting

    //  // Step 1. setting counts
    //  for(size_t i=0; i<n_particles; i++) {
    //    count[position_to_cube_id(scd, unsorted_particles[i])]++;
    //  }

    //  // Step 2. setting start index
    //  // this reuses the count array by summing up to a
    //  // continuous count
    // std::vector<SearchCube> retc;
    //  retc.reserve(scd.nt);

    //  size_t first = 0;
    //  size_t last = 0;
    //  for(size_t i=0; i<scd.nt; i++) {
    //    last = count[i];
    //    count[i] = first;

    //    retc.push_back({first, first+last});

    //    first += last;
    //  }

    SortedNeighbours ret {};
    //  ret.ids.reserve(unsorted_particles.size());

    // // Iterate search cubes
    //  const float maxDistanceSqr = scd.dx*scd.dx;
    //  std::vector<Point> retp(n_particles);
    //  NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
    //  SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

    //  size_t last_sid = 0;
    //  for(size_t uid=0; uid<n_particles; uid++) {
    //    // copy particle to new position in sorted array
    //    const size_t sid = position_to_cube_id(scd, unsorted_particles[uid]);
    //    size_t oid = count[sid]; // next
    //    retp[oid] = Point(unsorted_particles[uid]);
    //    const Point& opos = retp[oid];

    //    const auto ncids = lower_neighbour_cubes(sub, ncidsten, sid);

    //    // search for neighbours in lower search cube neighbours
    //    for (size_t ncid: ncids) {

    //      const size_t first_nc = retc[ncid].first;
    //      // search only up to count
    //      const size_t last_nc = count[ncid];

    //      // Step 3.3. set iterate neighbour cube particle
    //      for (size_t nid=first_nc; nid<last_nc; nid++) {

    //        const Point &npos = retp[nid];

    //        const float distanceSqr =
    //          squared_distance(opos, npos);

    //        if (distanceSqr < maxDistanceSqr) {
    //          ret.ids.push_back({oid, nid});
    //        };
    //      }
    //    }

    //    // Wenn sid sid < last sid particle backfill neighbour cube forward
    //    search if (sid < last_sid) {
    //      // TODO stimmen die uncids bei backfill ??
    //      const auto uncids = upper_neighbour_cubes(sub, ncidsten, sid);

    //      // search for neighbours in lower search cube neighbours
    //      for (size_t uncid: ncids) {

    //        const size_t first_unc = retc[uncid].first;
    //        // search only up to count
    //        const size_t last_unc = count[uncid];

    //        // Step 3.3. set iterate neighbour cube particle
    //        for (size_t nid=first_unc; nid<last_unc; nid++) {

    //          const Point &npos = retp[nid];

    //          const float distanceSqr =
    //            squared_distance(opos, npos);

    //          if (distanceSqr < maxDistanceSqr) {
    //            ret.ids.push_back({oid, nid});
    //          };
    //        }
    //      }
    //    }

    //    last_sid = sid;

    //    // Wenn sid voll dann owner cube search
    //    if (count[sid]==retc[sid].last) {
    //      owner_cube_search(retp, first, last, maxDistanceSqr, ret);
    //    };

    //    count[sid]++;

    //    // if oid was last particle in current sid search all neighbour
    //    // pairs in sid

    //    // find pairs in current sid
    //    // iterate over particles in sid with lower oid
    //    // size_t first = retc[sid].first;
    //    // for (size_t nid=oid; nid>first; nid--) {
    //    //   const Point &npos = retp[nid];

    //    //   const float distanceSqr =
    //    //     squared_distance(opos, npos);

    //    //   if (distanceSqr < maxDistanceSqr) {
    //    //     ret.ownId.push_back(oid);
    //    //     ret.neighId.push_back(nid);
    //    //   };
    //    // }
    //  }

    return ret;
}

// void owner_cube_search(first, last, ret,) {};
SortedNeighbours createNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> &pos,
    std::vector<SearchCube> &searchCubes) {}

SortedNeighbours createSTLNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> &pos,
    std::vector<SearchCube> &searchCubes,
    const std::vector<Facet_handle> &facets) {
    // Step 0 initialise return values

    SortedNeighbours ret {std::vector<NeighbourPair>(0),
                          std::vector<STLSurfaceDist>(0)};

    ret.ids.reserve(40 * pos.size());
    ret.dist.reserve(40 * pos.size());

#pragma omp parallel
    {
        const float maxDistanceSqr = scd.dx * scd.dx;

        NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
        SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

        // Step 1. get parent search cube
        SortedNeighbours ret_tmp {std::vector<NeighbourPair>(0),
                                  std::vector<STLSurfaceDist>(0)};

        ret_tmp.ids.reserve(40 * pos.size() / omp_get_num_threads());
        ret_tmp.dist.reserve(40 * pos.size() / omp_get_num_threads());
        // // serial
        // ret_tmp.ids.reserve(40 * pos.size());
        // ret_tmp.dist.reserve(40 * pos.size());

#pragma omp for nowait
        for (size_t sid = 0; sid < scd.nt; sid++) {

            // Step 2. test all particle in parent search cube
            // Step 2.1 set first particle of parent cube
            const size_t first = searchCubes[sid].first;
            const size_t last = searchCubes[sid].last;

            owner_cube_search(
                pos, first, last, maxDistanceSqr, facets, ret_tmp);

            const auto ncids = upper_neighbour_cubes(sub, ncidsten, sid);

            // Step 3. test all particles in neighbour search_cubes
            // Step 3.1. set neighbour search cube

            for (size_t ncid : ncids) {
                const size_t first_nc = searchCubes[ncid].first;
                const size_t last_nc = searchCubes[ncid].last;

                neighbour_cube_search(
                    pos,
                    first,
                    last,
                    first_nc,
                    last_nc,
                    maxDistanceSqr,
                    facets,
                    ret_tmp);
            }
        }

#pragma omp critical
        {
            ret.ids.insert(
                ret.ids.end(), ret_tmp.ids.begin(), ret_tmp.ids.end());
            ret.dist.insert(
                ret.dist.end(), ret_tmp.dist.begin(), ret_tmp.dist.end());
        }
    }

    return ret;
}
