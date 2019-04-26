#ifndef SEARCHCUBE2_H
#define SEARCHCUBE2_H

#include "SearchCubes.h"
#include <omp.h>
#include <math.h>
#include <random>

struct SearchCube2 {

    // bool empty;

    std::vector<size_t> neighCubes;

};

struct UnsortedSearchCubes {

    std::vector<bool> empty;

    std::vector<SearchCube2> search_cubes;

};

size_t position_zindex(
        float x, float y, float z,
        float min_x, float min_y, float min_z,
        size_t nx, size_t ny, size_t nz,
        float inv_dx, float dx
       ) {

    x = x - min_x;
    y = y - min_y;
    z = z - min_z;

    // const size_t ci{std::min(nx - 1, (size_t)(x * inv_dx))};
    // const size_t cj{std::min(ny - 1, (size_t)(y * inv_dx))};
    // const size_t ck{std::min(nz - 1, (size_t)(z * inv_dx))};
    const size_t ci{(size_t)(x * inv_dx)};
    const size_t cj{(size_t)(y * inv_dx)};
    const size_t ck{(size_t)(z * inv_dx)};
    const size_t cube_id{ci + nx * cj + nx * ny * ck};

    const float subx{x - ci * dx};  // fmod(x, dx);
    const float suby{y - cj * dx};  // fmod(y, dx);
    const float subz{z - ck * dx};  // fmod(z, dx);

    const size_t i{(size_t)(4.0 * subx * inv_dx)};
    const size_t j{(size_t)(4.0 * suby * inv_dx)};
    const size_t k{(size_t)(4.0 * subz * inv_dx)};

    // std::cout
    //     << "x " << p.x()
    //     << " y " << p.y()
    //     << " z " << p.z()
    //     << " sub x " << subx
    //     << " sub y " << suby
    //     << " sub z " << subz
    //     << " i " << i
    //     << " j " << j
    //     << " k " << k
    //     << " dx_" << dx_
    //     << "  cube_id " << cube_id
    //     << "  n_cubes_ " << n_cubes_[0]
    //     << std::endl;
    //
    const size_t sub{i + j * 4 + k * 16};
    return sub + 64*cube_id;
}

class SearchCubes2: public SearchCubesBase<SearchCube2> {

    // Nomenclature:
    //
    // pid - particle id, position in memory
    // nid - neighbour id, position in memory
    // pidu - particle id before sorting
    // sid  - search cube id in memory
    // snid - search cube neigbhour id in memory
    // ksid - kartesian search cube id
    //
    //
    // Some calculations
    //
    // on my notebook L1 = 32KB, 4B float = 8000
    //
    // fuer 512 Particle = 2K Pos
    //
    // Pro Partikel - 40 Nachbarn, 1/40 Suchwuerfel
    //
    // sids, kids    - nParticle = 8B*2*512  =   8K
    // first_, last_ - nCubes  = 8B*2*512/40 = 0.2K
    //
    // ParicleNeighbours
    // oridId, neighId - nNeighbours = 40*512*2*8B = 320K
    // In 32K passen ca 100 Partikel = 4^3 Nachbarschaften
    // 40*100*8
    //
    // 512 Particles Case
    //
    // Cacheline 64b - 16floats floats


    private:

        std::vector<size_t>         neighbourId_stencil_;

        bool keep_empty_cubes_;

        // move to ParticleNeighbours
        std::vector<Point>          sorted_points_;

        // maps sorted point id to corresponding sorted search cube id;
        std::vector<size_t>         sorted_points_search_cubes_;

        LeanParticleNeighbours      particle_neighbours_;

        // an array mapping kartesian natural ids to the
        // reduced set of search cubes
        // empty search cubes are indicated by its stop id
        std::vector<size_t>      sids_;

        std::vector<size_t>      kids_;

        size_t kartesian_cube_id_stop_;

        // for every search cube the first pid
        std::vector<size_t>      first_;

        // for every search cube the last pid
        std::vector<size_t>      last_;
        //
        // std::vector<size_t>      next;

    public:

        SearchCubes2(): SearchCubesBase() {};

        SearchCubes2(
                Logger logger,
                std::vector<Point>& points,
                const float dx,
                bool keep_empty_cubes=true,
                bool invoke_setup=true):
            SearchCubesBase(logger, points, dx),
            neighbourId_stencil_(compute_neighbourId_stencil()),
            keep_empty_cubes_(keep_empty_cubes)
    {
        logger_.info()
            << "Initialising the SearchCubes of "
            << n_cubes_[0] << "x" << n_cubes_[1] << "x" << n_cubes_[2]
            << "=" << tot_n_cubes_
            << " search cubes for "
            << n_points_ << " particles";
        //
        // kartesian id to signal a non existing search cube
        kartesian_cube_id_stop_ = tot_n_cubes_;

        // for every particle corresponding kartesian cube id
        // without physical locality
        //

        // initialize with kartesian_cube_id_stop_
        // sids_ = std::vector<size_t> (
        //         tot_n_cubes_,
        //         kartesian_cube_id_stop_
        //         );
        //
        // // initialize with kartesian_cube_id_stop_
        // kids_ = std::vector<size_t> (
        //         tot_n_cubes_,
        //         kartesian_cube_id_stop_
        //         );

        if (invoke_setup) {
            setup(points);
        }
    };

    void setup (
            std::vector<Point>& points
            ) {

        sorted_points_ = quickSortParticles(points);

        // size_t ctr = 0;
        // for ( auto p: sorted_points_) {
        //
        //     size_t cube_id = position_to_cube_id(p);
        //     size_t sub_id =  position_zindex(p);
        //
        //     std::cout
        //         << ctr
        //         << " c " << cube_id
        //         << " s " << sub_id
        //         << " x " << p.x()
        //         << " y " << p.y()
        //         << " z " << p.z()
        //         << std::endl;
        //
        //     ctr++;
        //
        // }
        //
        //
        // std::vector<std::vector<size_t>> pidu_kid (
        //         tot_n_cubes_
        //         );
        // //
        // //
        // // // get for every particle its corresponding search
        // // // cube, in particle order
        // //
        // setup_unsorted(points, pidu_kid);
        // //
        // // // So far the searchCubes_ are not fully initialised
        // // // since the .first and .last entries are not valid
        // // // hence particles are resorted in memory and .first and
        // // // .last are set
        // sorted_points_ = sortParticles(
        //     searchCubes_,
        //     pidu_kid,
        //     points,
        //     sorted_points_search_cubes_,
        //     sids_
        // );

        particle_neighbours_ = compute_particle_neighbours(
                sorted_points_,
                // sorted_points_search_cubes_,
                particle_neighbours_
                );
    };

    std::vector<Point> sortParticles(
        std::vector<SearchCube2>& searchCubes,
        const std::vector<std::vector<size_t>>& particles_in_cubes,
        const std::vector<Point>& unsorted_points,
        std::vector<size_t>& sorted_points_search_cubes,
        std::vector<size_t>& sorted_search_cubes_ids
    ) {
        logger_.info_begin()
            << "Restructuring particles";

        std::vector<Point> sorted_points;
        sorted_points.reserve(n_points_);
        sorted_points_search_cubes = std::vector<size_t> (n_points_, kartesian_cube_id_stop_);

        size_t particle_ctr=0;
        size_t searchCubeCtr = 0;
        first_ = std::vector<size_t> (tot_n_cubes_, n_points_);
        last_  = std::vector<size_t> (tot_n_cubes_, n_points_);
        for (auto& searchCube: searchCubes)
        {
                const auto& particles_in_cube = particles_in_cubes[searchCubeCtr];
                const size_t n_particles_in_cube = particles_in_cube.size();

                first_[searchCubeCtr] = particle_ctr;
                last_[searchCubeCtr]  = particle_ctr + n_particles_in_cube - 1;


                // move particles to sorted list
                for (size_t i=0; i < n_particles_in_cube; i++) {
                    // const size_t tmp_ptr_ctr = particle_ctr + i;
                    // logger_.verbose(99)
                    //     << " searchCubeCtr " k<< searchCubeCtr
                    //     << " particle_ctr " << tmp_ptr_ctr
                    //     << " particlePos " << point.x()
                    //     << " y " << point.y() << " z " << point.z();
                    sorted_points.emplace_back(unsorted_points[particles_in_cube[i]]);
                    // sorted_points_search_cubes[particle_ctr + i] = searchCubeCtr;
                }
                particle_ctr += n_particles_in_cube;
                searchCubeCtr++;
        }
        return sorted_points;
    };


    void setup_unsorted (
            const std::vector<Point> &points,
            std::vector<std::vector<size_t>>& pidu_kid
            ) {
            UnsortedSearchCubes tmp_searchCubes =
                construct_particle_search_cube_relation(
                    points,
                    pidu_kid
                );


            if (!keep_empty_cubes_) {
                logger_.info_begin()
                    << "Removing empty search cubes: ";

                for (size_t i=0; i<tot_n_cubes_; i++)
                {
                    if (tmp_searchCubes.empty[i] == false) {
                        searchCubes_.emplace_back(
                            std::move(tmp_searchCubes.search_cubes[i]));
                    }
                }
                const size_t tot_reduced_searchcubes = searchCubes_.size();

                logger_.info_end()
                    << "Remaining search cubes: " <<  tot_reduced_searchcubes;
            } else {
                searchCubes_ = tmp_searchCubes.search_cubes;
            }
    };


    // Getter

    const std::vector<size_t>& get_sorted_points_search_cubes() const {
        return sorted_points_search_cubes_;
    };

    const std::vector<Point>& get_sorted_points() const {
        return sorted_points_;
    };

    std::vector<Point>& get_sorted_points() { return sorted_points_; };

    const LeanParticleNeighbours& get_particle_neighbours() const {
        return particle_neighbours_;
    };

    LeanParticleNeighbours& get_particle_neighbours() {
        return particle_neighbours_;
    };

    // Setters

    void set_particles(std::vector<Point>& points) {
        sorted_points_ = points;
    };

    LeanParticleNeighbours& compute_particle_neighbours(
            const std::vector<Point>&  points,
            // const std::vector<size_t> searchcube_of_points,
            LeanParticleNeighbours&   particle_neighbours
            ) {
        logger_.info_begin() << "Computing particle pair distances ";

        // NOTE here space for the particle neighbours is reserved
        // otherwise the stack might be smashed
        // TODO dont hardcode the upper particle neighbour number
        logger_.info() << "Reserving memory";
        // particle_neighbours.nElems.reserve(points.size());
        particle_neighbours.nElems = std::vector<size_t> (points.size(), 0);
        particle_neighbours.neighId.reserve(points.size()*50);
        // std::fill(particle_neighbours.nElems.begin(),
        //           particle_neighbours.nElems.end(), 0);

        logger_.info() << "Looping";
        for (size_t pid=0; pid < n_points_; pid++) {


            set_particle_neighbours(
                    points,
                    pid,
                    dx_,
                    particle_neighbours);
        }

        size_t number_pairs = particle_neighbours.neighId.size();
        logger_.info_end() << "Found " << number_pairs << " particle pairs";

        return particle_neighbours;
    }

    UnsortedSearchCubes construct_particle_search_cube_relation(
            const std::vector<Point>& points,
            std::vector<std::vector<size_t>>& pidu_kid
            )
        {
            logger_.info_begin() << "Constructing particle search cube relationship";

            // Store if neighbours ids are already computed for given cube
            std::vector<bool> computed_neighbour_ids (
                    tot_n_cubes_,
                    false
                    );

            const size_t particle_stop_id = n_points_;

            // temporary vector of all search cubes
            // contiguous by natural kartesian id
            // containing empty search cubes
            std::vector<SearchCube2> searchCubes(
                tot_n_cubes_,
                {
                    std::vector<size_t> {}
                }
            );

            std::vector<bool> empty (tot_n_cubes_, true);

            for (size_t i=0; i<n_points_; i++) {

                const size_t cube_id = position_to_cube_id(points[i]);

                // NOTE cache misses if particles
                // and cube_id are not aligned
                pidu_kid[cube_id].push_back(i);

                // std::cout << "[DEBUG] Particle position " << i
                //     << " x " << points[i][0]
                //     << " y " << points[i][1]
                //     << " cube_id " << cube_id
                //     << std::endl;

                if (empty[cube_id]) {
                    // searchCubes[cube_id].id = cube_id;
                    empty[cube_id] = false;

                    // Set neighbour ids
                    if (!computed_neighbour_ids[cube_id]) {
                        for (const auto nid: compute_knid(cube_id)) {
                            // std::cout
                            //      << "[DEBUG] setting neighbour " << nid
                            //      << std::endl;
                            //
                            //   Cache write misses ?
                            searchCubes[nid].neighCubes.push_back(cube_id);
                        }
                        computed_neighbour_ids[cube_id] = true;
                    }

                    // std::cout << "[DEBUG] done nids " << std::endl;
                }
                // std::cout << "[DEBUG] done particle " << std::endl;
            }
            logger_.info_end();

            return {empty, searchCubes};
        }


    std::vector<Point>& quickSortParticles(
        std::vector<Point>& unsorted_points //,
    ) {
        // std::vector<Point> ret (unsorted_points);
        std::sort(unsorted_points.begin(), unsorted_points.end(),
                  [&](Point& a, Point& b) {
                      // TODO precompute z-index
                      const size_t id_a = position_zindex(
                              a.x(),
                              a.y(),
                              a.z(),
                              bound_box_.min().x(),
                              bound_box_.min().y(),
                              bound_box_.min().z(),
                              n_cubes_[0],
                              n_cubes_[1],
                              n_cubes_[2],
                              inv_dx_,
                              dx_)
                              ;
                      const size_t id_b = position_zindex(
                              b.x(),
                              b.y(),
                              b.z(),
                              bound_box_.min().x(),
                              bound_box_.min().y(),
                              bound_box_.min().z(),
                              n_cubes_[0],
                              n_cubes_[1],
                              n_cubes_[2],
                              inv_dx_,
                              dx_)
                              ;
                      // const size_t id_a = position_to_cube_id(a);
                      // const size_t id_b = position_to_cube_id(b);
                      return id_a < id_b;
                  });

        set_first_last(unsorted_points);
        return unsorted_points;
    };

    void set_first_last(std::vector<Point>& points) {
        first_ = std::vector<size_t> (tot_n_cubes_, n_points_);
        last_  = std::vector<size_t> (tot_n_cubes_, n_points_);

        for (size_t pid=0; pid < n_points_; pid++) {
            const size_t cube_id = position_to_cube_id(points[pid]);
            const size_t first = first_[cube_id];

            if (pid < first) {
                first_[cube_id] = pid;
                last_[cube_id]  = pid;
            }

            const size_t last  = last_[cube_id];
            if (pid > last) last_[cube_id] = pid;
        }

    };

    LeanParticleNeighbours& update_particle_neighbours(std::vector<Point>& sorted_points)
    {
        // TODO
        //  * Points must be resorted
        // #pragma omp parallel for
        if (particle_neighbours_.neighId.size() > 0) {
            particle_neighbours_.nElems.clear();
            particle_neighbours_.neighId.clear();
        }

         compute_particle_neighbours(
            sorted_points,
            particle_neighbours_
            );

         return particle_neighbours_;
    };

    const std::vector<size_t> compute_neighbourId_stencil () const {
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
        // search cube ids start from min(x) to max(z)
        // in a rhs coordinate system
        // - the right neighbour (x-axis) is +1
        // - the back neighbour  (y-axis) is +(n_cubes_[0])
        // - the upper neighbour (z-axis) is +(n_cubes_[0]*n_cubes[1])


        // TODO leave it a size_t, iterate only till 12, since
        // the stencil is symmetric anyway
        const size_t ny_nx {n_cubes_[0] * n_cubes_[1]};
        const size_t nx {n_cubes_[0]};

        return {
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

    }

    // convert kartesian cube id to individual coordinates
    const std::vector<size_t> id_to_i_j_k(const size_t id) const {
        const size_t n_cubes_x = n_cubes_[0];
        const size_t n_cubes_y = n_cubes_[1];
        const size_t n_cubes_z = n_cubes_[2];

        // TODO replace division here

        const size_t k{min(n_cubes_z - 1, id / n_cubes_x * n_cubes_y)};
        const size_t j{min(n_cubes_y - 1, id / n_cubes_x)};
        const size_t i{min(n_cubes_x - 1, id - j * n_cubes_x - k * n_cubes_x * n_cubes_y)};
        return {i, j, k};
    }

    // get non boundary neighbours from neighbour id, false - boundary, true - inner
    const std::vector<bool> masked_neighbourId_stencil(const size_t id) const {

        const size_t n_cubes_x = n_cubes_[0];
        const size_t n_cubes_y = n_cubes_[1];
        const size_t n_cubes_z = n_cubes_[2];

        const size_t k{min(n_cubes_z - 1, id / n_cubes_x * n_cubes_y)};
        const size_t j{min(n_cubes_y - 1, id / n_cubes_x)};
        const size_t i{min(n_cubes_x - 1, id - j * n_cubes_x - k * n_cubes_x * n_cubes_y)};

        std::vector<bool> ret (26, true);
        // TODO test if wrapping it in a loop causes performance degradation
        // TODO can probably generated with constexpr
        // if i == 0 mask left
        if  (i == 0) {
            for (int i : {0, 3, 6, 9, 12, 14, 17, 20, 23}) {
                ret[i] = false;
            }
        }
        // if j == 0 mask front
        if  (j == 0) {
            for (int i : {0, 1, 2, 9, 10, 11, 17, 18, 19}) {
                ret[i] = false;
            }
        }
        // if k == 0 mask bottom
        if  (k == 0) {
            for (int i : {0, 1, 2, 3, 4, 5, 6, 7, 8}) {
                ret[i] = false;
            }
        }

        if  (i == n_cubes_x - 1) {
            for (int i : {2, 5, 8, 11, 13, 16, 19, 22, 25}) {
                ret[i] = false;
            }
        }
        // if j == 0 mask back
        if  (j == n_cubes_y - 1) {
            for (int i : {6, 7, 8, 14, 15, 16, 23, 24, 25}) {
                ret[i] = false;
            }
        }
        // if k == 0 mask top
        if  (k == n_cubes_z - 1) {
            for (int i : {17, 18, 19, 20, 21, 22, 23, 24, 25}) {
                ret[i] = false;
            }
        }
        return ret;
    }

    const std::vector<size_t> compute_knid (
            const size_t ksid) const {

        // computes the kartesian neighbour ids for a given search cube
        // considers boundaries

        // TODO check for domain boundaries
        const std::vector<bool> neighbour_mask = masked_neighbourId_stencil(ksid);

        std::vector<size_t> neighbourIds;
        neighbourIds.reserve(26);

        // compute lower neighbour ids first
        for (size_t i = 0; i<neighbourId_stencil_.size(); i++)
        {
            if (neighbour_mask[12-i] && neighbourId_stencil_[i] <= ksid) {
                neighbourIds.emplace_back(ksid - neighbourId_stencil_[i]);
            }
        }

        // compute uppper neighbour ids
        for (size_t i = 0; i<neighbourId_stencil_.size(); i++)
        {
            if (neighbour_mask[i+13]) {
                neighbourIds.emplace_back(ksid + neighbourId_stencil_[i]);
            }
        }

        return neighbourIds;
    }

    const SearchCube2& get_searchCube_by_kartesianID (const size_t ksid) const {

        return  searchCubes_[sids_[ksid]];
    }

    const std::vector<SearchCube2>& get_searchCubes() const {
        return searchCubes_;
    };

    const std::vector<size_t>& get_first() const {
        return first_;
    };

    const std::vector<size_t>& get_last() const {
        return last_;
    };


    // Update Particle Neighbours as side effect
    // TODO pass search cube id first than get particles and search
    void  set_particle_neighbours(
            const std::vector<Point>& sorted_points,
            const size_t pid,
            const float maxDistance,
            LeanParticleNeighbours & particle_neighbours
        ) {
            const Point& point {sorted_points[pid]};
            const size_t sid = position_to_cube_id(point); //searchcube_of_points[oid];
            const float maxDistanceSquared = maxDistance*maxDistance;

            // TODO search all pairs in in parent cube first
            for (size_t nid = first_[sid]; nid <= last_[sid]; nid++) {

                // do nothing if neighbour is original particle
                if (pid == nid) continue;

                const float squared_distance_ = squared_distance(
                        point, sorted_points[nid]);

                if (squared_distance_ < maxDistanceSquared) {
                    particle_neighbours.nElems[pid]++;
                    particle_neighbours.neighId.push_back(nid);
                };
            }

            const size_t particle_stop_id = n_points_;
            // todo store knid
            for (const auto kncid: compute_knid(sid)) {

                // Hier wird ein mapping von zu snid gebraucht
                const size_t snid = kncid; // {kids_[kncid]};

                if (first_[snid] == n_points_) continue;

                for (size_t nid = first_[snid]; nid <= last_[snid]; nid++) {
                    const float squared_distance_ = squared_distance(
                            point, sorted_points[nid]);

                    if ( squared_distance_ < maxDistanceSquared) {

                        particle_neighbours.nElems[pid]++;
                        particle_neighbours.neighId.push_back(nid);
                    };
                }
            }
    }


    size_t position_to_cube_id(const Point& p) const {
        // TODO scale bounding box a bit
        // TODO test if speed up n_cubes are copied to a const size_t nx ...
        // std::cout << "[DEBUG] position_to_cube_id " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        const size_t i {std::min( (size_t) n_cubes_[0]-1, (size_t)((p.x() - bound_box_.min().x())*inv_dx_))};
        const size_t j {std::min( (size_t) n_cubes_[1]-1, (size_t)((p.y() - bound_box_.min().y())*inv_dx_))};
        const size_t k {std::min( (size_t) n_cubes_[2]-1, (size_t)((p.z() - bound_box_.min().z())*inv_dx_))};
        return i + n_cubes_[0]*j + n_cubes_[0]*n_cubes_[1]*k;
    }

    // size_t position_zindex(const Point p) const {
    //
    //     const float x = p.x() - bound_box_.min().x();
    //     const float y = p.y() - bound_box_.min().y();
    //     const float z = p.z() - bound_box_.min().z();
    //
    //     const size_t ci {std::min( (size_t) n_cubes_[0]-1, (size_t)(x*inv_dx_))};
    //     const size_t cj {std::min( (size_t) n_cubes_[1]-1, (size_t)(y*inv_dx_))};
    //     const size_t ck {std::min( (size_t) n_cubes_[2]-1, (size_t)(z*inv_dx_))};
    //     const size_t cube_id = ci + n_cubes_[0]*cj + n_cubes_[0]*n_cubes_[1]*ck;
    //
    //     const float particle_per_cube = 4.0;
    //
    //     const float subx = fmod(x, dx_);
    //     const float suby = fmod(y, dx_);
    //     const float subz = fmod(z, dx_);
    //
    //     const size_t i = 4.0 * subx*inv_dx_;
    //     const size_t j = 4.0 * suby*inv_dx_;
    //     const size_t k = 4.0 * subz*inv_dx_;
    //
    //     // std::cout
    //     //     << "x " << p.x()
    //     //     << " y " << p.y()
    //     //     << " z " << p.z()
    //     //     << " sub x " << subx
    //     //     << " sub y " << suby
    //     //     << " sub z " << subz
    //     //     << " i " << i
    //     //     << " j " << j
    //     //     << " k " << k
    //     //     << " dx_" << dx_
    //     //     << "  cube_id " << cube_id
    //     //     << "  n_cubes_ " << n_cubes_[0]
    //     //     << std::endl;
    //     //
    //     const size_t sub = i + j*4 + k*16;
    //     return sub + 64*cube_id;
    // }
};

#endif
