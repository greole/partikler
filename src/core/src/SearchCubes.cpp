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

#include "Logger.hpp" // for MSG, Logger
#include "Scalar.hpp" // for MSG, Logger
#include <algorithm>  // std::min, std::max
#include <math.h>     // for ceil

SearchCubeDomain
initSearchCubeDomain(const std::pair<Vec3, Vec3> &bound_box, Scalar dx) {

    // bounds are scaled a bit to avoid particles on exact boundary
    const Scalar domain_extrusion = 0.01 * dx; // in dx
    Scalar min_x = bound_box.first[0] - domain_extrusion;
    Scalar min_y = bound_box.first[1] - domain_extrusion;
    Scalar min_z = bound_box.first[2] - domain_extrusion;
    Scalar max_x = bound_box.second[0] + domain_extrusion;
    Scalar max_y = bound_box.second[1] + domain_extrusion;
    Scalar max_z = bound_box.second[2] + domain_extrusion;

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

    const Scalar idx = 1.0 / dx;
    return SearchCubeDomain {
        Vec3 {min_x, min_y, min_z},
        Vec3 {max_x, max_y, max_z},
        dx,
        idx,
        SubDivision {n_cubes_x, n_cubes_y, n_cubes_z},
        (size_t)n_cubes_x * (size_t)n_cubes_y * (size_t)n_cubes_z,
    };
}

std::pair<Vec3, Vec3> bounding_box(std::vector<Vec3> const &particles) {
    Scalar minx, miny, minz;
    Scalar maxx, maxy, maxz;
    minx = miny = minz = std::numeric_limits<Scalar>::max();
    maxx = maxy = maxz = std::numeric_limits<Scalar>::min();

    for (auto &el : particles) {

        minx = std::min(minx, el[0]);
        maxx = std::max(maxx, el[0]);

        miny = std::min(miny, el[1]);
        maxy = std::max(maxy, el[1]);

        minz = std::min(minz, el[2]);
        maxz = std::max(maxz, el[2]);
    }

    return std::pair<Vec3, Vec3> {{minx, miny, minz}, {maxx, maxy, maxz}};
}

// inplace functions
void owner_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const Scalar maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret) {

    for (size_t oid = first; oid < last; oid++) {
        Vec3 const &opos = pos[oid];

        // starts from first+oidctr since oidctr are already
        // tested against this nid (distance pairs)
        // +1 to avoid testing oid==nid
        for (size_t nid = oid + 1; nid < last; nid++) {
            Vec3 const &npos = pos[nid];

            if (opos == npos) {
                continue;
            }

            auto point_dist = npos - opos;
            Vec3 sd = {point_dist[0], point_dist[1], point_dist[2]};

            const Scalar distanceSqr = sd * sd;

            if (distanceSqr < maxDistanceSqr) {
                ret.push_back({{oid, nid}, sd});
            };
        }
    }
}

// std::array<bool, 27> vector_inner_owner_cube_search(
//     Point opos,
//     std::vector<Point> &pos,
//     size_t first,
//     size_t last,
//     Scalar maxDistanceSqr) {

//     std::array<bool, 27> mask {};

//     for (size_t i = 0; i < last - first; i++) {
//         mask[i] = false;
//     }

//     for (size_t nid = first; nid < last; nid++) {

//         // const Point npos = pos[nid];

//         const Scalar distanceSqr = squared_distance(opos, pos[nid]);

//         mask[nid - first] = distanceSqr < maxDistanceSqr;
//     }

//     return mask;
// }

void neighbour_cube_search(
    const std::vector<Vec3> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const Scalar maxDistanceSqr,
    std::vector<UnsortedNeighbour> &ret) {

    // Step 3.1. set pivot particle
    for (size_t oid = first; oid < last; oid++) {
        Vec3 const &opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for (size_t nid = first_nc; nid < last_nc; nid++) {
            Vec3 const &npos = pos[nid];

            auto point_dist = npos - opos;
            Vec3 sd = {point_dist[0], point_dist[1], point_dist[2]};

            const Scalar distanceSqr = sd * sd;

            if (distanceSqr < maxDistanceSqr) {
                ret.push_back({{oid, nid}, sd});
            };
        }
    }
}

SubDivision id_to_i_j_k(const size_t id, const SubDivision sub) {
    // TODO the nz-1 limiter seems unnecessary
    const size_t nx = sub.nx;
    const size_t ny = sub.ny;
    const size_t nz = sub.nz;
    const size_t nxny = nx * ny;
    const unsigned int k {(unsigned int)std::min(nz - 1, id / nxny)};
    const unsigned int j {(unsigned int)std::min(ny - 1, (id - k * nxny) / nx)};
    const unsigned int i {
        (unsigned int)std::min(nx - 1, id - j * nx - k * nxny)};
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
    const SearchCubeDomain scd, const std::vector<Vec3> &unsorted_particles) {

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

    std::vector<Vec3> retp(n_particles);
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

        // std::cout << " oid " << oid << " uid " << uid << " sid " << sid <<
        // std::endl;

        rets.push_back(oid);

        retp[oid] = Vec3 {unsorted_particles[uid]};
        count[sid]++;
    }

    return {retc, rets, retp};
}

SortedNeighbours createNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Vec3> &pos,
    std::vector<SearchCube> &searchCubes) {
    // Step 0 initialise return values

    SortedNeighbours ret {std::vector<NeighbourPair>(0), std::vector<Vec3>(0)};

    // own_id_min and size
    std::vector<std::pair<size_t, size_t>> order(omp_get_max_threads(), {0, 0});

    size_t tot_pairs {0};

#pragma omp parallel
    {
        const Scalar maxDistanceSqr = scd.dx * scd.dx;

        NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
        SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

        // Step 1. get parent search cube
        std::vector<UnsortedNeighbour> ret_tmp(0);
        // std::vector<NeighbourPair>(0),
        // std::vector<STLSurfaceDist>(0)};

        ret_tmp.reserve(40 * pos.size() / omp_get_num_threads());
        // ret_tmp.dist.reserve(40 * pos.size() / omp_get_num_threads());
        // // serial
        // ret_tmp.ids.reserve(40 * pos.size());
        // ret_tmp.dist.reserve(40 * pos.size());

#pragma omp for nowait
        for (size_t sid = 0; sid < scd.nt; sid++) {

            // Step 2. test all particle in parent search cube
            // Step 2.1 set first particle of parent cube
            const size_t first = searchCubes[sid].first;
            const size_t last = searchCubes[sid].last;

            owner_cube_search(pos, first, last, maxDistanceSqr, ret_tmp);

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
                    ret_tmp);
            }
        }

        // sort ret_tmp
        std::sort(
            ret_tmp.begin(),
            ret_tmp.end(),
            [](const auto &lhs, const auto &rhs) {
                return lhs.ids.ownId < rhs.ids.ownId;
            });

        // figure out the position at which the sorted ret_tmps are to be
        // inserted
        // 1. get total number of neighbour pairs
        // 2. push min to a vector and compute offsets
        // 3. resize the ret vector and insert according to offsets

        // after sorting this is the lowest owner id;
        int thread_num = omp_get_thread_num();

        size_t own_id_min = ret_tmp[0].ids.ownId;

        size_t particle_pairs_thread = ret_tmp.size();

        order[thread_num].first = own_id_min; // use this for start
        order[thread_num].second =
            particle_pairs_thread; // use this for offsets

        // #pragma omp critical
        //         {
        //             // for (auto& el: ret_tmp) {
        //             //     std::cout << omp_get_thread_num() << ": " <<
        //             el.ids.ownId << " "
        //             //               << el.ids.neighId << std::endl;
        //             // }
        //         }

        // wait till all done
#pragma omp barrier
#pragma omp single
        {
            // sort thread ids by lowest owner id
            std::sort(
                order.begin(), order.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });

            // compute start index for each thread
            for (auto &start : order) {
                tot_pairs += start.second;
                start.second = tot_pairs;
            }

            ret.ids.resize(tot_pairs);
            ret.dist.resize(tot_pairs);
        }

#pragma omp critical
        {
            // ret_tmp_item are presorted

            // get start
            size_t id {0};

            // find reordered thread id by lowest owner id
            for (size_t id_ = 0; id_ < order.size(); id_++) {
                if (order[id_].first == own_id_min) id = id_;
            }

            // copy elements from ret_tmp to final ret struct
            size_t start_idx = order[id].second - particle_pairs_thread;
            size_t end_idx = order[id].second;
            for (size_t idx = 0; idx < particle_pairs_thread; idx++) {
                size_t target_idx = idx + start_idx;
                ret.ids[target_idx] = ret_tmp[idx].ids;
                ret.dist[target_idx] = ret_tmp[idx].dist;
            }
        }
    }

    return ret;
}

std::ostream &operator<<(std::ostream &os, NeighbourPair const &n) {
    os << " own " << n.ownId << " neigh " << n.neighId;
    return os;
}
