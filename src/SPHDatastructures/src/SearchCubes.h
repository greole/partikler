#ifndef SEARCHCUBER_H
#define SEARCHCUBER_H

struct SearchCube {

    // kartesian id
    size_t id;

    bool empty;

    std::vector<size_t> neighCubes;

    // vector index of first particle
    size_t first;

    // vector index of last particle
    size_t last;
};

class SearchCubeTree {

    private:

        // domain geometry
        K::Iso_cuboid_3     bound_box_;
        float               dx_;
        std::vector<size_t> n_cubes_;
        size_t              tot_n_cubes_;

        // search cubes in kartesian order without empty search cubes
        // also contains the kartesian id of neighbouring search cube
        std::vector<SearchCube>         searchCubes_;

        std::vector<size_t> neighbourId_stencil_;

        Logger logger_;

        bool keep_empty_cubes_;

        std::vector<Point> sorted_points_;

        std::vector<size_t> sorted_points_search_cubes_;

        ParticleNeighbours particle_neighbours_;

        std::vector<size_t> sorted_search_cubes_ids_;

    public:

        SearchCubeTree():
            bound_box_(),
            dx_(0),
            n_cubes_({0, 0, 0}),
            tot_n_cubes_(n_cubes_[0] * n_cubes_[1] * n_cubes_[2])
        {};

        SearchCubeTree(
                const std::vector<Point>& points,
                const float dx,
                Logger logger,
                bool keep_empty_cubes=true):
            bound_box_(bounding_box(points.begin(), points.end()))  ,
            dx_(dx),
            n_cubes_({
                std::max(
                    (size_t) 1,
                    (size_t)ceil((bound_box_.max().x() - bound_box_.min().x())/dx)
                ),
                std::max(
                    (size_t) 1,
                    (size_t)ceil((bound_box_.max().y() - bound_box_.min().y())/dx)
                ),
                std::max(
                    (size_t) 1,
                    (size_t)ceil((bound_box_.max().z() - bound_box_.min().z())/dx)
                    )
                }),
            tot_n_cubes_(n_cubes_[0] * n_cubes_[1] * n_cubes_[2]),
            neighbourId_stencil_(compute_neighbourId_stencil()),
            logger_(logger),
            keep_empty_cubes_(keep_empty_cubes)
    {

        const size_t n_points {points.size()};

        logger_.info()
            << "Initialising the SearchCubeTree of "
            << n_cubes_[0] << "x" << n_cubes_[1] << "x" << n_cubes_[2]
            << "=" << tot_n_cubes_
            << " search cubes for "
            << n_points << " particles";
        //
        // kartesian id to signal a non existing search cube
        const size_t kartesian_cube_id_stop = tot_n_cubes_ + 1;

        // for every particle corresponding kartesian cube id
        // without physical locality
        std::vector<std::vector<size_t>> particles_kartesian_cube_id (
                tot_n_cubes_
                );

        // temporary vector of all search cubes
        std::vector<SearchCube> tmp_searchCubes(
            tot_n_cubes_,
            {
                kartesian_cube_id_stop,
                true,
                std::vector<size_t> {},
                0, 0
            }
        );

        std::vector<bool> computed_neighbour_ids (
                tot_n_cubes_,
                false
                );

        // TODO avoid copying
        sorted_search_cubes_ids_ = std::vector<size_t> (
                tot_n_cubes_,
                kartesian_cube_id_stop
                );

        // get for every particle its corresponding search
        // cube, in particle order
        construct_particle_search_cube_relation(
            n_points,
            points,
            particles_kartesian_cube_id,
            tmp_searchCubes,
            computed_neighbour_ids
        );

        removeEmptySearchCubes(
            n_points,
            tmp_searchCubes,
            particles_kartesian_cube_id,
            points,
            sorted_points_,
            sorted_points_search_cubes_,
            sorted_search_cubes_ids_,
            searchCubes_, // use = operator
            keep_empty_cubes_
        );

        // std::tuple<std::vector<size_t>, std::vector<size_t>, std::vector<float>> distances;
        compute_particle_distances(sorted_points_, particle_neighbours_);
    };


    // Getter

    const std::vector<size_t> &get_sorted_points_search_cubes() const {
        return sorted_points_search_cubes_;
    };

    const std::vector<SearchCube> &get_searchCubes() const {
        return searchCubes_;
    };

    const std::vector<Point> &get_sorted_points() const {
        return sorted_points_;
    };

    std::vector<Point> &get_sorted_points() {
        return sorted_points_;
    };

    const ParticleNeighbours &get_particle_distances() const {
        return particle_neighbours_;
    };

    ParticleNeighbours & get_particle_distances() {
        return particle_neighbours_;
    };

    void compute_particle_distances(
        const std::vector<Point> &sorted_points,
        ParticleNeighbours& particle_neighbours
        )
    {
        logger_.info_begin() << "Computing particle pair distances ";


        for (size_t pointID=0; pointID<sorted_points.size(); pointID++) {
            compute_particle_neighbours(
                    sorted_points,
                    pointID,
                    sorted_points_search_cubes_[pointID],
                    dx_,
                    particle_neighbours);
        }

        size_t number_pairs = particle_neighbours.origId.size();
        logger_.info_end() << "Found " << number_pairs << " particle pairs";

    };

    void construct_particle_search_cube_relation(
            const size_t n_points,
            const std::vector<Point> & points,
            std::vector<std::vector<size_t>>  & particles_kartesian_cube_id,
            std::vector<SearchCube> & tmp_searchCubes,
            std::vector<bool> &computed_neighbour_ids
            )
        {
            logger_.info_begin() << "Constructing particle search cube relationship";
            for (size_t i=0; i<n_points; i++) {
                const size_t cube_id = position_to_cube_id(points[i]);
                // std::cout << "[DEBUG] Particle position " << i
                //     << " x " << points[i][0]
                //     << " y " << points[i][1]
                //     << " cube_id " << cube_id
                //     << std::endl;

                particles_kartesian_cube_id[cube_id].push_back(i);

                if (tmp_searchCubes[cube_id].empty) {
                    tmp_searchCubes[cube_id].id = cube_id;
                    tmp_searchCubes[cube_id].empty = false;

                    // Set neighbour ids
                    if (!computed_neighbour_ids[cube_id]) {
                        for (const auto nid: computeNeighbourIds(cube_id)) {
                            // std::cout << "[DEBUG] setting neighbour " << nid << std::endl;
                            tmp_searchCubes[nid].neighCubes.push_back(cube_id);
                        }
                        computed_neighbour_ids[cube_id] = true;
                    }
                    // std::cout << "[DEBUG] done nids " << std::endl;
                }
                // std::cout << "[DEBUG] done particle " << std::endl;
            }
            logger_.info_end();
        }

    void removeEmptySearchCubes(
        const size_t n_points,
        const std::vector<SearchCube>& tmp_searchCubes,
        const std::vector<std::vector<size_t>>& particles_kartesian_cube_id,
        const std::vector<Point>& points,
        std::vector<Point>& sorted_points,
        std::vector<size_t>& sorted_points_search_cubes,
        std::vector<size_t>& sorted_search_cubes_ids,
        std::vector<SearchCube>& searchCubes,
        bool keep_empty_cubes
    ) {
        logger_.info_begin()
            << "Removing empty search cubes and restructuring particles";

        size_t particle_ctr=0;
        sorted_points.reserve(n_points);
        sorted_points_search_cubes.reserve(n_points);

        size_t searchCubeCtr = 0;
        for (auto tmp_searchCube: tmp_searchCubes)
        {

            if (tmp_searchCube.empty == false) {

                const size_t n_particles_in_cube =
                    particles_kartesian_cube_id[tmp_searchCube.id].size();

                tmp_searchCube.first = particle_ctr;
                tmp_searchCube.last  = particle_ctr+n_particles_in_cube-1;

                // move particles to sorted list
                for (size_t i=0; i < n_particles_in_cube; i++) {
                    Point point {points[particles_kartesian_cube_id[tmp_searchCube.id][i]]};
                    // const size_t tmp_ptr_ctr = particle_ctr + i;
                    // logger_.verbose(99)
                    //     << " searchCubeCtr " << searchCubeCtr
                    //     << " particle_ctr " << tmp_ptr_ctr
                    //     << " particlePos " << point.x() << " y " << point.y() << " z " << point.z();
                    sorted_points.push_back(point);
                    sorted_points_search_cubes[particle_ctr + i] = searchCubeCtr;
                }
                particle_ctr += n_particles_in_cube;
                sorted_search_cubes_ids[tmp_searchCube.id] = searchCubeCtr;
                searchCubes.push_back(std::move(tmp_searchCube));
                searchCubeCtr++;
            }
            // if (keep_empty_cubes && tmp_searchCube.empty == true)  {
            //     sorted_search_cubes_ids[tmp_searchCube.id] = searchCubeCtr;
            //     searchCubes.push_back(std::move(tmp_searchCube));
            //     searchCubeCtr++;
            // }
        }
        const size_t tot_reduced_searchcubes = searchCubes.size();

        logger_.info_end()
            << "Remaining search cubes: " <<  tot_reduced_searchcubes;
    };

    ParticleNeighbours& update_particle_neighbours(std::vector<Point>& sorted_points) {
        // TODO
        //  * Points must be resorted
        // #pragma omp parallel for
        if (particle_neighbours_.origId.size() > 0) {
            particle_neighbours_.origId.clear();
            particle_neighbours_.neighId.clear();
            particle_neighbours_.distances.clear();
            particle_neighbours_.normalised_distances.clear();
            particle_neighbours_.squared_length.clear();
        }

         compute_particle_distances(
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

    const SearchCube get_searchCube_by_kartesianID (const size_t kartesian_cube_id) const {

        auto it = std::find_if(
            searchCubes_.begin(),
            searchCubes_.end(),
            [&kartesian_cube_id](const SearchCube& searchCube) {
                return kartesian_cube_id==searchCube.id;}
        );

        return *it;
    }

    void  compute_particle_neighbours(
            const std::vector<Point>& sorted_points,
            const size_t pointID,
            const size_t searchCubeId,
            const float maxDistance,
            ParticleNeighbours & particle_neighbours
        ) {
            const Point& point {sorted_points[pointID]};

            const SearchCube& parent_search_cube {searchCubes_[searchCubeId]};

            // NOTE here space for the particle neighbours is reserved
            // otherwise the stack might be smashed
            // TODO dont hardcode the upper particle neighbour number
            //
            particle_neighbours.origId.reserve(sorted_points.size()*40);
            particle_neighbours.neighId.reserve(sorted_points.size()*40);
            particle_neighbours.distances.reserve(sorted_points.size()*40);
            particle_neighbours.normalised_distances.reserve(sorted_points.size()*40);
            particle_neighbours.squared_length.reserve(sorted_points.size()*40);

            // std::cout << "[DEBUG] search parent cube " << std::endl;
            // TODO search all pairs in in parent cube first


            for (
                size_t particle_id = parent_search_cube.first;
                particle_id <= parent_search_cube.last;
                particle_id++) {


                // if (!bloom[neighId])     continue;
                if (pointID == particle_id) continue;

                const float squared_distance_ = squared_distance(
                        point, sorted_points[particle_id]);

                // std::cout
                //     << "[DEBUG] "
                //     << "point1 " << p[0] << " " << p[1] << " " << p[2]
                //     << "point2 " << p2[0] << " " << p2[1] << " " << p2[2]
                //     << " maxDistance*maxDistance " << maxDistance*maxDistance
                //     << " distance " << distance
                //     << std::endl;
                //

                if (squared_distance_ < maxDistance*maxDistance) {
                    // TODO use move semantics here, eg. emplace back
                    // this probably needs a move operator for the CGAL Vector
                    const float distance = sqrt(squared_distance_);

                    const Vector distance_vect {
                       sorted_points[particle_id].x() - point.x(),
                       sorted_points[particle_id].y() - point.y(),
                       sorted_points[particle_id].z() - point.z()
                    };

                    const Vector normalised_distance_vect {
                        distance_vect.x()/distance,
                        distance_vect.y()/distance,
                        distance_vect.z()/distance
                    };

                    // std::cout << "[DEBUG] push_back" << std::endl;
                    particle_neighbours.origId.push_back(pointID);
                    particle_neighbours.neighId.push_back(particle_id);
                    particle_neighbours.distances.push_back(distance_vect);
                    particle_neighbours.normalised_distances.push_back(normalised_distance_vect);
                    particle_neighbours.squared_length.push_back(distance);
                };
            }

            const std::vector<size_t> search_cube_neighbours =
                parent_search_cube.neighCubes;

            // std::cout << "[DEBUG] search neighbour cubes for sid "
            //           << searchCubeId << std::endl;
            for (const auto nid: search_cube_neighbours) {
                // std::cout << "[DEBUG] nid " << nid << std::endl;
                const size_t sorted_nid {sorted_search_cubes_ids_[nid]};
                // std::cout << "[DEBUG] sorted_nid " << sorted_nid << std::endl;
                for (
                    size_t particle_id = searchCubes_[sorted_nid].first;
                    particle_id <= searchCubes_[sorted_nid].last;
                    particle_id++) {


                    const float squared_distance_ = squared_distance(
                            point, sorted_points[particle_id]);

                    if ( squared_distance_ < maxDistance*maxDistance) {
                        const float distance = sqrt(squared_distance_);
                        const Vector distance_vect {
                             sorted_points[particle_id].x() - point.x(),
                             sorted_points[particle_id].y() - point.y(),
                             sorted_points[particle_id].z() - point.z()
                        };

                        const Vector normalised_distance_vect {
                            distance_vect.x()/distance,
                            distance_vect.y()/distance,
                            distance_vect.z()/distance
                        };

                        particle_neighbours.origId.push_back(pointID);
                        particle_neighbours.neighId.push_back(particle_id);
                        particle_neighbours.distances.push_back(distance_vect);
                        particle_neighbours.normalised_distances.push_back(normalised_distance_vect);
                        particle_neighbours.squared_length.push_back(distance);
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
