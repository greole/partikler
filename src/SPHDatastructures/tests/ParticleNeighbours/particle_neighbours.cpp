#include "CGALTYPEDEFS.h"
#include "gtest/gtest.h"
#include "SPHDatastructures.h"
#include <vector>

class SimpleParticleNeighbour :
    public ::testing::Test
{

    private:

        ParticleNeighbours particle_neighbours_;

    public:

        ParticleNeighbours& get_particle_neighbours() {return particle_neighbours_;};

        // Setup a 4 particle neighbour relation
        // 0 - {0., 0., 0.}  1 - {1., 0., 0.} 2 - {0., 1., 0.} 3 - {0., 0., 1.}
        SimpleParticleNeighbour( ) {
            particle_neighbours_ = ParticleNeighbours {
                {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3},
                {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2},
                {
                    {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},     // 0
                    {-1., 0., 0.}, {-1., 1., 0.}, {-1., 0., 1.},  // 1
                    {0., -1., 0.}, {1., -1., 0.}, {0., -1., 1.},  // 2
                    {0., 0., -1.}, {1., 0., -1.}, {0., 1., -1.}  // 3
                },
                { // TODO correct these
                    {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},     // 0
                    {-1., 0., 0.}, {-1., 1., 0.}, {-1., 0., 1.},  // 1
                    {0., -1., 0.}, {1., -1., 0.}, {0., -1., 1.},  // 2
                    {0., 0., -1.}, {1., 0., -1.}, {0., 1., -1.}  // 3
                },
                {1., 2., 3., 0., 2., 3., 0., 1., 3., 0., 1., 2.}
            };
        }

        ~SimpleParticleNeighbour( ) {
        }

};

TEST_F (SimpleParticleNeighbour, add_ab0) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values
    std::vector<float> src {0.0, 0.0, 0.0, 0.0};
    std::vector<float> dst; // = std::vector(12, 0.0);

    dst = particle_neighbours.add(src);

    ASSERT_EQ (dst.size(), 12);

    for (auto d: dst) {
        ASSERT_EQ (d,  0.0);
    }
}

TEST_F (SimpleParticleNeighbour, add_ab1) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values
    std::vector<float> src {1.0, 1.0, 1.0, 1.0};
    std::vector<float> dst; // = std::vector(12, 0.0);

    dst = particle_neighbours.add(src);

    ASSERT_EQ (dst.size(), 12);

    for (auto d: dst) {
        ASSERT_EQ (d,  2.0);
    }
}

TEST_F (SimpleParticleNeighbour, add_ab2) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length // Make this a separate test
    //  * correct values
    std::vector<float> src {1.0, 1.0, 0.0, 0.0};
    std::vector<float> dst;

    dst = particle_neighbours.add(src);

    ASSERT_EQ (dst.size(), 12);

    std::vector<float> trgt {
        2.0, 1.0, 1.0,
        2.0, 1.0, 1.0,
        1.0, 1.0, 0.0,
        1.0, 1.0, 0.0
    };

    for (size_t ctr=0; ctr<dst.size(); ctr++) {
        ASSERT_EQ (dst[ctr],  trgt[ctr]);
    }
}


TEST_F (SimpleParticleNeighbour, sub_ab0) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values
    std::vector<float> src {0.0, 0.0, 0.0, 0.0};
    std::vector<float> dst; // = std::vector(12, 0.0);

    dst = particle_neighbours.subtract(src);

    ASSERT_EQ (dst.size(), 12);

    for (auto d: dst) {
        ASSERT_EQ (d,  0.0);
    }
}

TEST_F (SimpleParticleNeighbour, sub_ab1) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values
    std::vector<float> src {1.0, 1.0, 1.0, 1.0};
    std::vector<float> dst; // = std::vector(12, 0.0);

    dst = particle_neighbours.subtract(src);

    ASSERT_EQ (dst.size(), 12);

    for (auto d: dst) {
        ASSERT_EQ (d,  0.0);
    }
}

TEST_F (SimpleParticleNeighbour, getb1) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values
    std::vector<float> src {1.0, 2.0, 3.0, 4.0};
    std::vector<float> dst; // = std::vector(12, 0.0);

    dst = particle_neighbours.get_b(src);

    ASSERT_EQ (dst.size(), 12);

    for (auto d: dst) {
        std::cout << d << std::endl;
        // ASSERT_EQ (d,  0.0);
    }
}

TEST_F (SimpleParticleNeighbour, sum) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values

    std::vector<float> src {
        1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0
    };

    std::vector<float> dst (4, 0.0);

    particle_neighbours.sum(src, dst);

    ASSERT_EQ (dst.size(), 4);

    for (auto d: dst) {
        ASSERT_EQ (d,  3.0);
    }
}

TEST_F (SimpleParticleNeighbour, sum2) {

    ParticleNeighbours& particle_neighbours = get_particle_neighbours();

    // TODO tests:
    //  * correct length
    //  * correct values

    std::vector<float> src {
        1.0, 2.0, 3.0,
        1.0, 1.0, 2.0,
        2.0, 1.0, 2.0,
        3.0, 2.0, 1.0
    };

    std::vector<float> dst (4, 0.0);

    particle_neighbours.sum(src, dst);

    ASSERT_EQ (dst.size(), 4);

    std::vector<float> trgt {
        6.0, 4.0, 5.0, 6.0
    };

    for (size_t ctr=0; ctr<dst.size(); ctr++) {
        ASSERT_EQ (dst[ctr],  trgt[ctr]) << ctr;
    }

    // for (auto d: dst) {
    //     ASSERT_EQ (d,  6.0);
    // }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
