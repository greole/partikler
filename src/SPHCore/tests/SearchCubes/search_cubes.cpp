#include "CGAL_TYPEDEFS.h"
#include "gtest/gtest.h"
#include "SearchCubes.h"
#include "computes.h"
#include "io.h"
#include <vector>

// TODO check if global setup and tear down
// allows to an instance
class KartesianPointGrid :
    public ::testing::Test
{

    public:

        // Setup a 10 by 10 particle grid
        // 4 by 4 particles in a cube
        KartesianPointGrid( ) {
            for(int j = 0; j<10; j++) {
                for(int i = 0; i<10; i++) {
                    points.push_back(Point((float) i, (float) j, 0.0));
                }
            }
            search_cubes = new SearchCubeTree(points, 3.33, RunTime(4));
        }

        ~KartesianPointGrid( ) {
            delete search_cubes;
        }

        // void SetUp( ) {
        // }

        std::vector<Point> points;

        SearchCubeTree *search_cubes;
};


TEST_F(KartesianPointGrid, KartesianIDtoIJKconversion) {
    size_t ctr = 0;
    for (size_t j: {0, 1, 2}) {
        for (size_t i: {0, 1, 2}) {
            const std::vector<size_t> stencil = search_cubes->id_to_i_j_k(ctr);
            const std::vector<size_t> target = {i, j, 0};
            ASSERT_EQ (target, stencil ) << ctr;
            ctr++;
        }
    }
}

TEST_F(KartesianPointGrid, DefaultStencil) {
    // test if the default stencil is correct
    const std::vector<size_t> stencil = search_cubes->compute_neighbourId_stencil();

    const std::vector<size_t> target_stencil {
          1,   2,   3,   4,  5,  6,  7,  8,  9, 10, 11, 12, 13};

    for (size_t i = 0, e = target_stencil.size(); i != e; i++) {
        ASSERT_EQ (stencil[i],  target_stencil[i]) << i;
    }
}

TEST_F(KartesianPointGrid, ParticleInSearchCubes) {

    std::vector<size_t> n_particles_searchcube
    {
        16, 12, 12,
        12, 9, 9,
        12, 9, 9
    };

    size_t particle_ctr = 0;
    size_t cube_ctr = 0;

    for (auto np: n_particles_searchcube) {
        for (int i=0; i<np; i++) {
                size_t searchCubeId = search_cubes->get_sorted_points_search_cubes()[particle_ctr];
                ASSERT_EQ (cube_ctr, searchCubeId)
                    << " particle_ctr: " << particle_ctr
                    << " cube_ctr: " << cube_ctr ;
                particle_ctr++;
        }
        cube_ctr++;
    }
}

TEST_F(KartesianPointGrid, ParticleNeighbours) {

    const ParticleNeighbours &particle_neighbours
        {search_cubes->get_particle_distances()};

    const std::vector<Point> sorted_points {search_cubes->get_sorted_points()};

    size_t ctr = 0;
    for (auto & src_particle_id: particle_neighbours.origId) {
        const size_t neigbour_particle_id {particle_neighbours.neighId[ctr]};
        const float distance {particle_neighbours.squared_length[ctr]};
        // std::cout
        //     << "src_particle_id: " << src_particle_id
        //     << " neigbour_particle_id " << neigbour_particle_id
        //     << " distance " << distance
        //     << std::endl;
        ctr++;
    }
}

TEST_F(KartesianPointGrid, SearchCubeNeighbours) {

    const std::vector<SearchCube> &sorted_search_cubes = search_cubes->get_searchCubes();

    ASSERT_EQ (sorted_search_cubes.size(),  9);

    // number of particles in each search cube
    std::vector<size_t> n_particles {16, 12, 12, 12, 9, 9, 12, 9, 9};

    // neighbour ids
    std::vector<std::vector<size_t>> nids = {
         {1, 3, 4}, {0, 2, 3, 4, 5}, {1, 4, 5},
         {0, 1, 4, 6, 7}, {0, 1, 2, 3, 5, 6, 7, 8}, {1, 2, 4, 7, 8},
         {3, 4, 7}, {3, 4, 5, 6, 8}, {4, 5, 7}
    };

    size_t particle_ctr = 0;
    for (int sid = 0; sid < 9; sid++) {

        // all search cubes should be populated
        ASSERT_EQ (sorted_search_cubes[sid].empty,  false);

        // test if number of neihbours is correct
        ASSERT_EQ (sorted_search_cubes[sid].neighCubes.size(),  nids[sid].size())
            << " sid " << sid;

        ASSERT_EQ (sorted_search_cubes[sid].first, particle_ctr);
        particle_ctr += n_particles[sid];
        ASSERT_EQ (sorted_search_cubes[sid].last,  particle_ctr - 1);

        size_t nid_ctr = 0;
        for (auto nid: nids[sid]) {
            ASSERT_EQ (sorted_search_cubes[sid].neighCubes[nid_ctr],  nid);
            nid_ctr++;
        }
    }
}

// test if the boundaries are correctly masked
TEST_F(KartesianPointGrid, StencilMaskBottomFrontLeft) {

    // TODO test all boundaries
    const size_t id = 0;
    const std::vector<bool> mask = search_cubes->masked_neighbourId_stencil(id);

    // complete top and bottom are masked
    std::vector<bool> target (26, false);

    target[13] = true;
    target[15] = true;
    target[16] = true;

    for (size_t i = 0, e = mask.size(); i != e; i++) {
        ASSERT_EQ (mask[i],  target[i]) << i;
    }
}

// test if the boundaries are correctly masked
TEST_F(KartesianPointGrid, StencilMaskTopBackRight) {

    // TODO test all boundaries
    const size_t id = 26;
    const std::vector<bool> mask = search_cubes->masked_neighbourId_stencil(id);

    std::vector<bool> target (26, false); // complete bottom is masked

    target[9] = true;
    target[10] = true;
    target[12] = true;

    for (size_t i = 0, e = mask.size(); i != e; i++) {
        ASSERT_EQ (mask[i],  target[i]) << i;
    }
}

TEST_F(KartesianPointGrid, Computes) {

    Kernel kernel;
    const auto particle_neighbours {search_cubes->get_particle_distances()};
    kernel.W    = std::move(std::vector<float>  (particle_neighbours.origId.size(), 0));
    kernel.dWdx = std::move(std::vector<Vector> (particle_neighbours.origId.size(), {0, 0, 0}));

    RunTime runTime(99);
    compute_kernel(runTime, 2.0, particle_neighbours, kernel);

    const std::vector<Point> points = search_cubes->get_sorted_points();
    const size_t n_points = points.size();
    std::vector<float> densities (n_points, 0.0);

    compute_density(runTime, particle_neighbours, kernel, densities);

    for(auto dens: densities) {
        std::cout << "[DEBUG] dens "  << dens << std::endl;
    }

    // writeData_SPH("grid_test", 0, points, densities);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
