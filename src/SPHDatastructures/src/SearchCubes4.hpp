#ifndef SEARCHCUBE4_H
#define SEARCHCUBE4_H

struct SearchCube4 {

    // kartesian id
    size_t id;

    bool empty;

    std::vector<size_t> neighCubes;

    // vector index of first particle
    size_t first;

    // vector index of last particle
    size_t last;
};


class SearchCubes4: public SearchCubesBase<SearchCube4> {


    // Naive Implementation
    // 2D:
    // =====
    // Number neighbourhood particles nbp = 40
    // Number particles per Searchcube = npb*A_cube/A_circ = nbp * 4/Pi = 1.3*npb
    // Number of particles per Mega search cube = nbp*1.3*9 = 11.7
    // 3D
    // =====
    // Number particles per Searchcube = npb*V_cube/V_sphere = nbp * 6/Pi = 1.9*npb
    // Number of particles per Mega search cube = nbp*1.9*27 = 51.3
    //
    // ULLR
    // 2D:
    // Misses = 4/9*9*1.3 = 5.2
    // 3D:
    // Misses = 8/27*27*1.9 = 15.2

    private:

        std::vector<size_t>         neighbourId_stencil_;

        bool keep_empty_cubes_;

        // move to ParticleNeighbours
        std::vector<Point>          sorted_points_;

        // maps sorted point id to corresponding sorted search cube id;
        std::vector<size_t>         sorted_points_search_cubes_;

        LeanParticleNeighbours      particle_neighbours_;

        std::vector<size_t>         sorted_search_cubes_ids_;

        size_t kartesian_cube_id_stop_;

    public:

        SearchCubes4(): SearchCubesBase() {};

        SearchCubes4(
                Logger logger,
                const std::vector<Point>& points,
                const float dx,
                bool keep_empty_cubes=true):
            SearchCubesBase(logger, points, dx),
            neighbourId_stencil_(compute_neighbourId_stencil()),
            keep_empty_cubes_(keep_empty_cubes)
    {

        const size_t n_points {points.size()};

        logger_.info()
            << "Initialising the SearchCubes of "
            << n_cubes_[0] << "x" << n_cubes_[1] << "x" << n_cubes_[2]
            << "=" << tot_n_cubes_
            << " search cubes for "
            << n_points << " particles";
        //
        // kartesian id to signal a non existing search cube
        kartesian_cube_id_stop_ = tot_n_cubes_ + 1;

        // for every particle corresponding kartesian cube id
        // without physical locality
        std::vector<std::vector<size_t>> particles_kartesian_cube_id (
                tot_n_cubes_
                );


        // TODO avoid copying
        // an array mapping kartesian natural ids to the
        // reduced set of search cubes
        // empty search cubes are indicated by its stop id
        sorted_search_cubes_ids_ = std::vector<size_t> (
                tot_n_cubes_,
                kartesian_cube_id_stop_
                );

        /* std::cout << "sizeof" << sizeof(particles_kartesian_cube_id)  << std::endl; */
        /* std::cout << "alignof" << alignof(particles_kartesian_cube_id)  << std::endl; */

        // get for every particle its corresponding search
        // cube, in particle order
        {
            std::vector<SearchCube4> tmp_searchCubes =
                construct_particle_search_cube_relation(
                    points,
                    particles_kartesian_cube_id
                );

            if (!keep_empty_cubes) {
                logger_.info_begin()
                    << "Removing empty search cubes: ";

                size_t cubeCtr = 0;
                size_t nonEmptyCubeCtr = 0;

                for (auto tmp_searchCube: tmp_searchCubes)
                {
                    if (tmp_searchCube.empty == false) {
                        searchCubes_.push_back(tmp_searchCube);
                        sorted_search_cubes_ids_[cubeCtr] = nonEmptyCubeCtr;
                        nonEmptyCubeCtr++;
                    }
                    cubeCtr++;
                }
                const size_t tot_reduced_searchcubes = searchCubes_.size();

                logger_.info_end()
                    << "Remaining search cubes: " <<  tot_reduced_searchcubes;
            } else {
                searchCubes_ = tmp_searchCubes;
            }
        }

        // So far the searchCubes_ are not fully initialised
        // since the .first and .last entries are not valid
        // hence particles are resorted in memory and .first and
        // .last are set
        sorted_points_ = sortParticles(
            searchCubes_,
            particles_kartesian_cube_id,
            points,
            sorted_points_search_cubes_,
            sorted_search_cubes_ids_
        );

        // size_t ctr = 0;
        // for ( auto p: sorted_points_) {
        //
        //
        //     std::cout
        //         << ctr
        //         // << " c " << cube_id
        //         // << " s " << sub_id
        //         << " x " << p.x()
        //         << " y " << p.y()
        //         << " z " << p.z()
        //         << std::endl;
        //
        //     ctr++;
        //
        // }

        // std::tuple<std::vector<size_t>, std::vector<size_t>, std::vector<float>> distances;
        // compute_particle_distances(sorted_points_, particle_neighbours_);
        compute_particle_neighbours(
                points,
                sorted_points_search_cubes_,
                particle_neighbours_
                );
    };

    // Getter

    const std::vector<size_t> &get_sorted_points_search_cubes() const {
        return sorted_points_search_cubes_;
    };

    const std::vector<Point> &get_sorted_points() const {
        return sorted_points_;
    };

    std::vector<Point> &get_sorted_points() {
        return sorted_points_;
    };

    const LeanParticleNeighbours &get_particle_distances() const {
        return particle_neighbours_;
    };

    LeanParticleNeighbours & get_particle_distances() {
        return particle_neighbours_;
    };

    LeanParticleNeighbours & get_particle_neighbours() {
        return particle_neighbours_;
    };


    LeanParticleNeighbours compute_particle_neighbours(
            const std::vector<Point>  points,
            const std::vector<size_t> searchcube_of_points,
            LeanParticleNeighbours&       particle_neighbours
            ) {
        logger_.info_begin() << "Computing particle pair distances ";

        // NOTE here space for the particle neighbours is reserved
        // otherwise the stack might be smashed
        // TODO dont hardcode the upper particle neighbour number
        logger_.info() << "Reserving memory";
        particle_neighbours.nElems = std::vector<size_t> (points.size(), 0);
        particle_neighbours.neighId.reserve(points.size()*50);

        logger_.info() << "Looping";
        for (size_t oid=0; oid < n_points_; oid++) {

            const size_t searchcube_of_point = searchcube_of_points[oid];

            set_particle_neighbours(
                    points,
                    oid,
                    searchcube_of_point,
                    dx_,
                    particle_neighbours);
        }

        size_t number_pairs = particle_neighbours.neighId.size();
        logger_.info_end() << "Found " << number_pairs << " particle pairs";

        return particle_neighbours;
    }

    std::vector<SearchCube4> construct_particle_search_cube_relation(
            const std::vector<Point>& points,
            std::vector<std::vector<size_t>>& particles_kartesian_cube_id
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
            std::vector<SearchCube4> searchCubes(
                tot_n_cubes_,
                {
                    kartesian_cube_id_stop_,
                    true,
                    std::vector<size_t> {},
                    particle_stop_id, 0
                }
            );


            for (size_t i=0; i<n_points_; i++) {

                const size_t cube_id = position_to_cube_id(points[i]);
                particles_kartesian_cube_id[cube_id].push_back(i);

                // std::cout << "[DEBUG] Particle position " << i
                //     << " x " << points[i][0]
                //     << " y " << points[i][1]
                //     << " cube_id " << cube_id
                //     << std::endl;

                if (searchCubes[cube_id].empty) {
                    searchCubes[cube_id].id = cube_id;
                    searchCubes[cube_id].empty = false;


                    // if (i < tmp_searchCubes[cube_id].first) {
                    //     tmp_searchCubes[cube_id].first = i;
                    // }
                    //
                    // if (i > tmp_searchCubes[cube_id].last) {
                    //     tmp_searchCubes[cube_id].last = i;
                    // }

                    // Set neighbour ids
                    if (!computed_neighbour_ids[cube_id]) {
                        for (const auto nid: computeNeighbourIds(cube_id)) {
                            // std::cout
                            //      << "[DEBUG] setting neighbour " << nid
                            //      << std::endl;
                            searchCubes[nid].neighCubes.push_back(cube_id);
                        }
                        computed_neighbour_ids[cube_id] = true;
                    }
                    // std::cout << "[DEBUG] done nids " << std::endl;
                }
                // std::cout << "[DEBUG] done particle " << std::endl;
            }
            logger_.info_end();

            return searchCubes;
        }

    std::vector<Point> sortParticles(
        std::vector<SearchCube4>& searchCubes,
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

        size_t particle_ctr = 0;
        size_t searchCubeCtr = 0;
        for (auto& searchCube: searchCubes)
        {
                const auto& particles_in_cube = particles_in_cubes[searchCube.id];
                const size_t n_particles_in_cube = particles_in_cube.size();

                searchCube.first = particle_ctr;
                searchCube.last  = particle_ctr + n_particles_in_cube - 1;

                // move particles to sorted list
                for (size_t i=0; i < n_particles_in_cube; i++) {
                    // const size_t tmp_ptr_ctr = particle_ctr + i;
                    // logger_.verbose(99)
                    //     << " searchCubeCtr " << searchCubeCtr
                    //     << " particle_ctr " << tmp_ptr_ctr
                    //     << " particlePos " << point.x()
                    //     << " y " << point.y() << " z " << point.z();
                    sorted_points.emplace_back(unsorted_points[particles_in_cube[i]]);
                    sorted_points_search_cubes[particle_ctr + i] = searchCubeCtr;
                }
                particle_ctr += n_particles_in_cube;
                searchCubeCtr++;
        }
        return sorted_points;
    };
    //
    // ParticleNeighbours& update_particle_neighbours(std::vector<Point>& sorted_points) {
    //     // TODO
    //     //  * Points must be resorted
    //     // #pragma omp parallel for
    //     if (particle_neighbours_.neighId.size() > 0) {
    //         particle_neighbours_.origId.clear();
    //         particle_neighbours_.neighId.clear();
    //         particle_neighbours_.distances.clear();
    //         particle_neighbours_.normalised_distances.clear();
    //         particle_neighbours_.squared_length.clear();
    //     }
    //
    //      compute_particle_neighbours(
    //         sorted_points,
    //         sorted_points_search_cubes_,
    //         particle_neighbours_
    //         );
    //
    //      return particle_neighbours_;
    // };

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
    const std::vector<size_t> id_to_i_j_k (const size_t id) const {
        // id =  i + n_cubes_[0]*j + n_cubes_[0]*n_cubes_[1]*k;
        // k = id // n_cubes_[0]*n_cubes_[1]
        const size_t k {min( n_cubes_[2] - 1, id / n_cubes_[0]*n_cubes_[1])};
        const size_t j {min( n_cubes_[1] - 1, id / n_cubes_[0])};
        const size_t i {min( n_cubes_[0] - 1, id - j*n_cubes_[0] - k*n_cubes_[0]*n_cubes_[1])};
        // std::cout
        //     << "[DEBUG] id " << id
        //     << " n_cubes_[0] " << n_cubes_[0]
        //     << " n_cubes_[1] " << n_cubes_[1]
        //     << " n_cubes_[2] " << n_cubes_[2]
        //     << " k  " << k
        //     << " j  " << j
        //     << " i  " << i
        //     << std::endl;
        return {i, j, k};
    }

    // get non boundary neighbours from neighbour id, false - boundary, true - inner
    const std::vector<bool> masked_neighbourId_stencil(const size_t id) const {
        const std::vector<size_t> ijk = id_to_i_j_k(id);
        std::vector<bool> ret (26, true);
        // TODO test if wrapping it in a loop causes performance degradation
        // TODO can probably generated with constexpr
        // if i == 0 mask left
        if  (ijk[0] == 0) {
            for (int i : {0, 3, 6, 9, 12, 14, 17, 20, 23}) {
                ret[i] = false;
            }
        }
        // if j == 0 mask front
        if  (ijk[1] == 0) {
            for (int i : {0, 1, 2, 9, 10, 11, 17, 18, 19}) {
                ret[i] = false;
            }
        }
        // if k == 0 mask bottom
        if  (ijk[2] == 0) {
            for (int i : {0, 1, 2, 3, 4, 5, 6, 7, 8}) {
                ret[i] = false;
            }
        }

        if  (ijk[0] == n_cubes_[0] - 1) {
            for (int i : {2, 5, 8, 11, 13, 16, 19, 22, 25}) {
                ret[i] = false;
            }
        }
        // if j == 0 mask back
        if  (ijk[1] == n_cubes_[1] - 1) {
            for (int i : {6, 7, 8, 14, 15, 16, 23, 24, 25}) {
                ret[i] = false;
            }
        }
        // if k == 0 mask top
        if  (ijk[2] == n_cubes_[2] - 1) {
            for (int i : {17, 18, 19, 20, 21, 22, 23, 24, 25}) {
                ret[i] = false;
            }
        }
        return ret;
    }

    const std::vector<size_t> computeNeighbourIds (
            const size_t id) const {

        // TODO check for domain boundaries
        // std::cout << " [DEBUG] computing lower neighbours for cube " << id << std::endl;

        const std::vector<bool> neighbour_mask = masked_neighbourId_stencil(id);

        std::vector<size_t> neighbourIds;

        // compute lower neighbour ids first
        for (size_t i = 0; i<neighbourId_stencil_.size(); i++)
        {
            if (neighbour_mask[12-i]) {
                // TODO remove this check once neighbour_mask works properly
                // check if nid underflows size_t
                // std::cout << "[DEBUG] testing lower nid id " << id << " i " << i  << std::endl;
                if (neighbourId_stencil_[i] <= id) {
                    const size_t nid {id - neighbourId_stencil_[i]};
                    // std::cout << "[DEBUG] pushing lower nid" << nid << std::endl;
                    // TODO REMOVE this check, why is it necessary?
                    neighbourIds.push_back(nid);
                }
            }
        }

        // std::cout << " [DEBUG] computing upper neighbours for cube " << id << std::endl;
        // compute lower neighbour ids first
        for (size_t i = 0; i<neighbourId_stencil_.size(); i++)
        {
            if (neighbour_mask[i+13]) {
                const size_t nid {neighbourId_stencil_[i] + id};
                // std::cout << "[DEBUG] nid" << nid << std::endl;
                // std::cout << "[DEBUG] pushing upper nid" << nid << std::endl;
                // TODO REMOVE this check, why is it necessary?
                neighbourIds.push_back(nid);
            }
        }
        return neighbourIds;
    }

    const SearchCube4 get_searchCube_by_kartesianID (const size_t kartesian_cube_id) const {

        auto it = std::find_if(
            searchCubes_.begin(),
            searchCubes_.end(),
            [&kartesian_cube_id](const SearchCube4& searchCube) {
                return kartesian_cube_id==searchCube.id;}
        );

        return *it;
    }

    // Update Particle Neighbours as side effect
    void  set_particle_neighbours(
            const std::vector<Point>& sorted_points,
            const size_t oid,
            const size_t searchCubeId,
            const float maxDistance,
            LeanParticleNeighbours& particle_neighbours
        ) {
            const Point& point {sorted_points[oid]};

            const SearchCube4& parent_search_cube {searchCubes_[searchCubeId]};

            for (
                size_t nid = parent_search_cube.first;
                nid <= parent_search_cube.last;
                nid++) {

                // do nothing if neighbour is original particle
                if (oid == nid) continue;

                const float squared_distance_ = squared_distance(
                        point, sorted_points[nid]);

                if (squared_distance_ < maxDistance*maxDistance) {
                    particle_neighbours.nElems[oid]++;
                    particle_neighbours.neighId.push_back(nid);
                };
            }

            const std::vector<size_t> search_cube_neighbours =
                parent_search_cube.neighCubes;

            for (const auto nid: search_cube_neighbours) {
            /* for (const auto nid: computeNeighbourIds(searchCubeId)) { */
                const size_t sorted_nid {sorted_search_cubes_ids_[nid]};
                for (
                    size_t nid = searchCubes_[sorted_nid].first;
                    nid <= searchCubes_[sorted_nid].last;
                    nid++) {

                    const float squared_distance_ = squared_distance(
                            point, sorted_points[nid]);

                    if ( squared_distance_ < maxDistance*maxDistance) {
                        particle_neighbours.nElems[oid]++;
                        particle_neighbours.neighId.push_back(nid);
                    };
                }
            }
    }

    const size_t position_to_cube_id(const Point& p) const {
        // TODO scale bounding box a bit
        // TODO test if speed up n_cubes are copied to a const size_t nx ...
        // std::cout << "[DEBUG] position_to_cube_id " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        const size_t i {std::min( (size_t) n_cubes_[0]-1, (size_t)((p.x() - bound_box_.min().x())/dx_))};
        const size_t j {std::min( (size_t) n_cubes_[1]-1, (size_t)((p.y() - bound_box_.min().y())/dx_))};
        const size_t k {std::min( (size_t) n_cubes_[2]-1, (size_t)((p.z() - bound_box_.min().z())/dx_))};
        return i + n_cubes_[0]*j + n_cubes_[0]*n_cubes_[1]*k;
    }
};


#endif
