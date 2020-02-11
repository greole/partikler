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

#include "CGALSearchCubes.hpp"

#include "Logger.hpp" // for MSG, Logger
#include "Scalar.hpp" // for MSG, Logger

// inplace functions
void stl_owner_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const Scalar maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret) {

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

            const Scalar distanceSqr = sd.len * sd.len;

            if (distanceSqr < maxDistanceSqr) {
                ret.push_back({{oid, nid}, sd});
            };
        }
    }
}

void stl_neighbour_cube_search(
    const std::vector<Point> &pos,
    const size_t first,
    const size_t last,
    const size_t first_nc,
    const size_t last_nc,
    const Scalar maxDistanceSqr,
    const std::vector<Facet_handle> &facets,
    std::vector<STLUnsortedNeighbour> &ret) {

    // Step 3.1. set pivot particle
    for (size_t oid = first; oid < last; oid++) {
        const Point &opos = pos[oid];

        // Step 3.3. set iterate neighbour cube particle
        for (size_t nid = first_nc; nid < last_nc; nid++) {
            const Point &npos = pos[nid];

            STLSurfaceDist sd =
                compute_STLSurface_dist(opos, npos, facets[oid], facets[nid]);

            const Scalar distanceSqr = sd.len * sd.len;

            if (distanceSqr < maxDistanceSqr) {
                ret.push_back({{oid, nid}, sd});
            };
        }
    }
}

STLSortedNeighbours createSTLNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> &pos,
    std::vector<SearchCube> &searchCubes,
    const std::vector<Facet_handle> &facets) {
    // Step 0 initialise return values

    STLSortedNeighbours ret {std::vector<NeighbourPair>(0),
                             std::vector<STLSurfaceDist>(0)};

    // own_id_min and size
    std::vector<std::pair<size_t, size_t>> order(omp_get_max_threads(), {0, 0});

    size_t tot_pairs {0};

#pragma omp parallel
    {
        const Scalar maxDistanceSqr = scd.dx * scd.dx;

        NeighbourIdHalfStencil ncidsten(scd.n.nx, scd.n.ny);
        SubDivision sub {scd.n.nx, scd.n.ny, scd.n.nz};

        // Step 1. get parent search cube
        std::vector<STLUnsortedNeighbour> ret_tmp(0);
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

            stl_owner_cube_search(
                pos, first, last, maxDistanceSqr, facets, ret_tmp);

            const auto ncids = upper_neighbour_cubes(sub, ncidsten, sid);

            // Step 3. test all particles in neighbour search_cubes
            // Step 3.1. set neighbour search cube

            for (size_t ncid : ncids) {
                const size_t first_nc = searchCubes[ncid].first;
                const size_t last_nc = searchCubes[ncid].last;

                stl_neighbour_cube_search(
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

std::pair<Vec3, Vec3> cgal_bounding_box(std::vector<Point> const &particles) {
    Scalar minx, miny, minz;
    Scalar maxx, maxy, maxz;
    minx, miny, minz = std::numeric_limits<Scalar>::max();
    maxx, maxy, maxz = std::numeric_limits<Scalar>::min();

    for (auto &el : particles) {

        minx = min(minx, (Scalar)el.x());
        maxx = max(maxx, (Scalar)el.x());

        miny = min(miny, (Scalar)el.y());
        maxy = max(maxy, (Scalar)el.y());

        minz = min(minz, (Scalar)el.z());
        maxz = max(maxz, (Scalar)el.z());
    }

    return std::pair<Vec3, Vec3> {{minx, miny, minz}, {maxx, maxy, maxz}};
}
