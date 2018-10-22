
struct SearchCube {

    // kartesian id
    size_t id;

    bool empty;

    std::vector<size_t> neighCubes;

    // vector index of first particle
    size_t first;

    // vector index of last particle
    size_t last;

    // std::vector<Point>       points;
};

class SearchCubeTree {

    private:

        // domain geometry
        const K::Iso_cuboid_3           bound_box_;
        const float                     dx_;
        const std::vector<size_t>       n_cubes_;
        const size_t                    tot_n_cubes_;


        // search cubes in kartesian order without empty search cubes
        // also contains the kartesian id of neighbouring search cube
        std::vector<SearchCube>         searchCubes;

        std::vector<size_t> neighbourId_stencil_;

        std::vector<Point> sorted_points;

    public:

        SearchCubeTree(const std::vector<Point>& points, const float dx):
            bound_box_(bounding_box(points.begin(), points.end())),
            dx_(dx),
            n_cubes_({
               (size_t)ceil((bound_box_.max().x() - bound_box_.min().x())/dx),
               (size_t)ceil((bound_box_.max().y() - bound_box_.min().y())/dx),
               (size_t)ceil((bound_box_.max().z() - bound_box_.min().z())/dx)
                }),
            tot_n_cubes_(n_cubes_[0] * n_cubes_[1] * n_cubes_[2])
    {

        std::cout
            << "Initialising the SearchCubeTree for "
            << n_cubes_[0] << "x" << n_cubes_[1] << "x" << n_cubes_[2]
            << "=" << tot_n_cubes_
            <<  " search cubes" << std::endl;

        const size_t n_points {points.size()};

        // kartesian id to signal a non existing search cube
        const size_t kartesian_cube_id_stop = tot_n_cubes_ + 1;

        // for every particle corresponding kartesian cube id
        // without physical locality
        // std::vector<size_t> particles_kartesian_cube_id (
        //         n_points, kartesian_cube_id_stop
        //         );
        //
        std::vector<std::vector<size_t>> particles_kartesian_cube_id (
                tot_n_cubes_
                );

        SearchCube empty = {
                    kartesian_cube_id_stop,
                    true,
                    std::vector<size_t> {},
                    0, 0
        };

        std::vector<SearchCube>  tmp_searchCubes(
                tot_n_cubes_, empty
                );

        // get for every particle its corresponding search
        // cube, in particle order
        for (size_t i=0; i<n_points; i++) {
            const size_t cube_id = position_to_cube_id(points[i]);
            particles_kartesian_cube_id[i].push_back(cube_id);

            tmp_searchCubes[cube_id].empty = false;

            // Set neighbour ids
            for (const auto nid: computeNeighbourIds(cube_id)) {
                tmp_searchCubes[nid].neighCubes.push_back(cube_id);
            }

        }

        std::cout
            << "Removing empty search cubes \n"
            << "and restructuring particles"
            << std::endl;

        size_t particle_ctr=0;

        for (auto & tmp_searchCube: tmp_searchCubes) {
            if (tmp_searchCube.empty == false) {

                const size_t n_particles_in_cube =
                    particles_kartesian_cube_id[tmp_searchCube.id].size();

                tmp_searchCube.first = particle_ctr;
                tmp_searchCube.last  = particle_ctr+n_particles_in_cube;
                particle_ctr += n_particles_in_cube;

                // move particles to sorted list
                for (size_t i=0; i < n_particles_in_cube; i++) {
                    Point point {points[particles_kartesian_cube_id[tmp_searchCube.id][i]]};
                    sorted_points.push_back(point);
                }

                searchCubes.push_back(std::move(tmp_searchCube));
            }
        }

        std::cout << searchCubes.size() << std::endl;

        std::tuple<std::vector<size_t>, std::vector<size_t>, std::vector<float>> distances;
        for (const auto& point: sorted_points) {
            compute_particle_neighbours(point, dx_, distances);
        }

        std::cout << std::get<2>(distances)[0] << std::endl;

    };

    const std::vector<size_t> neighbourId_stencil () const {

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
        //
        // TODO Implement Search Cube Boundaries
        // On the domain boundaries:
        //
        //  Bottom        Middle        Top
        //  ^ y           ^ y           ^ y
        //  |05 06 --     |11 12 --     |17 18 -- // back
        //  |03 04 --     |09 10 --     |15 16 -- // middle
        //  |01 02 --     |07 08 --     |13 14 -- // front
        //  |-------> x   |-------> x   |-------> x
        //  left  right                           n_cubes_[0]=3
        //
        //  -- stop value
        // - the right neighbour (x-axis) is --
        // - the back neighbour  (y-axis) is +(n_cubes_[0])-1
        // - the upper neighbour (z-axis) is +(n_cubes_[0]*n_cubes[1])-1

        std::vector<size_t> offsets (26);

        const size_t ny_nx {n_cubes_[0] * n_cubes_[1]};
        const size_t nx {n_cubes_[0]};

        size_t ctr = 0;
        size_t stencil_id;
        for (int k: {-1, 0, +1}) {
            for (int j: {-1, 0, +1}) {
                for (int i: {-1, 0, +1}) {
                    if ( ctr != 13) { // do nothing for the centre cube
                        // after the centre cube the stencil index has to be
                        // adapted
                        (ctr > 12) ? stencil_id = ctr -1 : stencil_id = ctr;

                        offsets[stencil_id] = k*ny_nx + j*nx + i;
                        ctr++;
                    }
                }
            }
        }

        return offsets;
    }

    const std::vector<size_t> computeNeighbourIds (
            const size_t id) const {

        std::vector<size_t> neighbourIds;
        for (
            auto it = neighbourId_stencil_.begin();
            it != neighbourId_stencil_.end();
            it++)
        {
            neighbourIds.push_back(*it + id);
        }
        return neighbourIds;
    }

    const SearchCube get_searchCube_by_kartesianID (const size_t kartesian_cube_id) const {

        auto it = std::find_if(
            searchCubes.begin(),
            searchCubes.end(),
            [&kartesian_cube_id](const SearchCube& searchCube) {
                return kartesian_cube_id==searchCube.id;}
        );

        return *it;
    }

    void  compute_particle_neighbours(
            const Point& point, const float maxDistance,
            std::tuple<std::vector<size_t>, std::vector<size_t>, std::vector<float>> distances
            ) {
            // TODO probably better to pass point id;

            std::pair<std::vector<size_t>, std::vector<float>>  ret;

            const size_t search_cube_kartesian_id = position_to_cube_id(point);

            const SearchCube parent_search_cube =
                get_searchCube_by_kartesianID(search_cube_kartesian_id);


            for (
                size_t particle_id = parent_search_cube.first;
                particle_id < parent_search_cube.last;
                particle_id++) {

                const float distance = squared_distance(
                        point, sorted_points[particle_id]);

                if ( distance < maxDistance*maxDistance) {
                    std::get<0>(distances).push_back(particle_id);
                    std::get<1>(distances).push_back(particle_id);
                    std::get<2>(distances).push_back(distance);
                };

            }

            const std::vector<size_t> search_cube_neighbours =
                parent_search_cube.neighCubes;

            for (const auto nid: search_cube_neighbours) {
                for (
                    size_t particle_id = parent_search_cube.first;
                    particle_id < parent_search_cube.last;
                    particle_id++) {

                    const float distance = squared_distance(
                            point, sorted_points[particle_id]);

                    if ( distance < maxDistance*maxDistance) {
                        std::get<0>(distances).push_back(particle_id);
                        std::get<1>(distances).push_back(particle_id);
                        std::get<2>(distances).push_back(distance);
                    };
                }
            }
    }

    const size_t position_to_cube_id(const Point& p) const {
        const size_t i {(size_t)((p.x() - bound_box_.min().x())/dx_)};
        const size_t j {(size_t)((p.y() - bound_box_.min().y())/dx_)};
        const size_t k {(size_t)((p.z() - bound_box_.min().z())/dx_)};
        return i + n_cubes_[0]*j + n_cubes_[0]*n_cubes_[1]*k;
    }
};

