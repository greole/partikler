#include "SearchCubes.hpp"
#include "Helper.hpp"
#include "gtest/gtest.h"

#include <vector>

bool found_particle_pair(
    size_t i, size_t j, size_t k, std::vector<UnsortedNeighbour> &ret) {
    bool found_i = (ret[k].ids.ownId == i);
    bool found_j = (ret[k].ids.neighId == j);
    return (found_i and found_j);
}

// TODO needs facets
TEST(OwnerCubeSearch, FindsNeighbours) {
    // number of particles in each direction
    const size_t n_particles = 3;
    std::vector<Vec3> points =
        create_uniform_particle_cube(n_particles); // evtl andere  funktion

    const size_t tot_n_particles = n_particles * n_particles * n_particles;

    std::vector<UnsortedNeighbour> ret {};
    ret.reserve(tot_n_particles * points.size() / 2);

    owner_cube_search(points, 0, tot_n_particles, 5.0, ret);

    // Should find at least some neighbours
    ASSERT_GT(ret.size(), 0);

    // Neighbourhood should be found once at max
    ASSERT_LT(ret.size(), tot_n_particles * tot_n_particles);

    // When each particle has all particles as neighbour
    // the number of particle pairs is n-1+n-2+...
    ASSERT_EQ(ret.size(), (tot_n_particles - 1) / 2 * tot_n_particles);

    // Every particle needs to be found either as owner or neighbour
    // TODO Refactor
    // for (size_t i=0; i<tot_n_particles; i++) {
    //   bool present =
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end()
    //     ||
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end();

    //   ASSERT_TRUE(present) << i;
    // }

    // Each combination of (i,j) or (j,i) needs to present

    for (size_t i = 0; i < tot_n_particles; i++) {
        size_t present = 0;
        for (size_t j = i + 1; j < tot_n_particles; j++) {
            for (size_t k = 0; k < ret.size(); k++) {

                if (found_particle_pair(i, j, k, ret) or
                    found_particle_pair(j, i, k, ret)) {
                    present++;
                }
            }
            ASSERT_GT(present, 0)
                << "Particle pair " << i << " " << j << " missing";
            ASSERT_LT(present, 2) << "Found particle pair " << i << " " << j
                                  << "  more than once";
            present = 0;
        }
    }
}

TEST(NeighbourCubeSearch, FindsNeighbours) {
    // number of particles in each direction
    const size_t n_particles = 3;
    std::vector<Vec3> points = create_uniform_particle_cube(n_particles);
    std::vector<Vec3> neighbour_points =
        create_uniform_particle_cube(n_particles);

    // shift neighbour particles x-position by 1
    for (size_t i = 0; i < neighbour_points.size(); i++) {
        auto point = neighbour_points[i];
        neighbour_points[i] =
            Vec3 {(Scalar)(point[0] + 1.0), point[1], point[2]};
    }

    // append neighbour points to owner points
    for (auto point : neighbour_points) {
        points.push_back(point);
    }

    const size_t tot_n_particles = n_particles * n_particles * n_particles;

    std::vector<UnsortedNeighbour> ret {};
    ret.reserve(tot_n_particles * points.size() / 2);

    neighbour_cube_search(
        points,
        0,
        tot_n_particles,
        tot_n_particles,
        2 * tot_n_particles,
        10.0,
        ret);

    // Should find at least some neighbours
    ASSERT_GT(ret.size(), 0);

    // When each particle has all particles as neighbour
    // the number of particle pairs is n_particles**2
    ASSERT_EQ(ret.size(), tot_n_particles * tot_n_particles);

    // Every particle needs to be found either as owner or neighbour
    // TODO refactor
    // for (size_t i=0; i<tot_n_particles; i++) {
    //   bool present =
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end()
    //     ||
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end();

    //   ASSERT_TRUE(present) << i;
    // }

    // every particle of the owner cube should find all neighbour cube
    // particles as neighbours
    for (size_t i = 0; i < tot_n_particles; i++) {
        size_t present = 0;
        for (size_t j = tot_n_particles; j < 2 * tot_n_particles; j++) {
            for (size_t k = 0; k < ret.size(); k++) {

                if (found_particle_pair(i, j, k, ret) or
                    found_particle_pair(j, i, k, ret)) {
                    present++;
                }
            }
            ASSERT_GT(present, 0)
                << "Particle pair " << i << " " << j << " missing";
            ASSERT_LT(present, 2) << "Found particle pair " << i << " " << j
                                  << "  more than once";
            present = 0;
        }
    }
}

TEST(ParticleNeighbourSearch, createNeighbours) {
    // Same tests as for owner_cube_search apply

    // number of particles in each direction
    const size_t n_particles = 12;
    std::vector<Vec3> points = create_uniform_particle_cube(n_particles);

    const size_t tot_n_particles = n_particles * n_particles * n_particles;

    const size_t n_cubes = 64;
    const size_t tot_n_particles_search_cube = tot_n_particles / 64;
    float cube_dx = 1.0 / ((float)n_particles);

    auto bb = bounding_box(points);
    const SearchCubeDomain scd = initSearchCubeDomain(bb, 0.25);

    SortedParticles sp = countingSortParticles(scd, points);

    // TODO searchCubes are created by counting sort
    std::vector<SearchCube> searchCubes(1, {0, points.size()});

    SortedNeighbours ret = createNeighbours(scd, sp.particles, sp.searchCubes);

    // Should find at least some neighbours
    ASSERT_GT(ret.ids.size(), 0);

    // Neighbourhood should be found once at max
    ASSERT_LT(ret.ids.size(), tot_n_particles * tot_n_particles);

    // Neighbourhood should be found once at max are
    // tot_n_particles_search_cube**2*n_cubes
    // ASSERT_LT(ret.ownId.size(), tot_n_particles*tot_n_particles/64);

    // When each particle has all particles as neighbour
    // the number of particle pairs are approx n-1+n-2+...

    const size_t estimate_n_neighbours = tot_n_particles * 40;

    ASSERT_GT(ret.ids.size(), estimate_n_neighbours * 3 / 4);
    ASSERT_LT(ret.ids.size(), estimate_n_neighbours * 4 / 3);

    // Every particle needs to be found either as owner or neighbour
    // TODO Refactor
    // for (size_t i=0; i<tot_n_particles; i++) {
    //   bool present =
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end()
    //     ||
    //     std::find(ret.ids.begin(), ret.ids.end(), i) != ret.ids.end();

    //   ASSERT_TRUE(present) << i;
    // }
}

struct SearchCubesTestData {
    float cube_dx;
    std::vector<Vec3> points;
    size_t expected_n_cubes;
    std::vector<size_t> expected_neighbour_ids;

    SearchCubesTestData(
        float cube_dx,
        std::vector<Vec3> points,
        size_t expected_n_cubes,
        std::vector<size_t> expected_neighbour_ids)
        : cube_dx(cube_dx), points(points), expected_n_cubes(expected_n_cubes),
          expected_neighbour_ids(expected_neighbour_ids) {};
};

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
