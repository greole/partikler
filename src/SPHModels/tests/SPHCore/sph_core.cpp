#include "CGALTYPEDEFS.h"
#include "gtest/gtest.h"
#include "SPHDatastructures.h"
#include "SPHio.h"
#include "SPHCore.h"
#include <vector>

class KartesianPointGrid :
    public ::testing::Test
{

    public:

        // Setup a 10 by 10 particle grid
        // 4 by 4 particles in a0 cube
        const size_t n = 40;
        KartesianPointGrid( ) {
            for(size_t j = 0; j<n; j++) {
                for(size_t i = 0; i<n; i++) {
                    float x = 0;
                    float y = 0;
                    x = ((float)(rand() % n))/(3.*(float) n);
                    y = ((float)(rand() % n))/(3.*(float) n);
                    x += (float) i;
                    y += (float) j;
                    if (( i > 18 && i < 22 ) && ( j > 18 && j < 22 )) { continue; }
                    points.push_back(Point(x, y, 0.0));
                }
            }
            search_cubes = new SearchCubeTree(points, 3.33, Logger(0));
        }

        ~KartesianPointGrid( ) {
            delete search_cubes;
        }

        // void SetUp( ) {
        // }

        std::vector<Point> points;

        SearchCubeTree *search_cubes;
};

void set_boundary(std::vector<int>& ids, std::vector<Point>& pos) {

    for(size_t ctr = 0; ctr<pos.size(); ctr++) {

        if (
            pos[ctr].x() < 5. || pos[ctr].x() > 34. ||
            pos[ctr].y() < 5. || pos[ctr].y() > 34.
            ) {
            ids[ctr] = 0;
        }

    }
}

TEST_F(KartesianPointGrid, Computes) {

    Logger logger {99};
    RunTime runTime {logger, false};

    const float dx = 1.0;
    runTime.set_dict("dx", dx);
    runTime.set_dict("n_timesteps", (size_t) 5000);

    std::vector<Point> & points = search_cubes->get_sorted_points();

    SPHPointField positions = runTime.set_particle_positions(points);

    ParticleNeighbours& particle_neighbours = runTime.initialize_particle_neighbours();

    std::vector<int>    & type = runTime.create_uniform_int_field(1.0,  "ID").get_field();
    std::vector<int>    & sid  = runTime.create_uniform_int_field(1.0, "SID").get_field();
    std::vector<float>  & rho = runTime.create_uniform_field(1.0      ,"rho").get_field();
    std::vector<float>  &   p = runTime.create_uniform_field(1.01e5   ,  "p").get_field();
    std::vector<Vector> & dp  = runTime.create_uniform_field({0, 0, 0}, "dp").get_field();
    std::vector<Vector> & dnu = runTime.create_uniform_field({0, 0, 0},"dnu").get_field();
    std::vector<Vector> &  du = runTime.create_uniform_field({0, 0, 0}, "du").get_field();
    std::vector<Vector> &   u = runTime.create_uniform_field({0, 0, 0},  "u").get_field();

    set_boundary(type, points);

    const float dt = 5e-05;
    for (;!runTime.end(); runTime++) {
        std::cout << "Timestep " << runTime.get_timeStep() << std::endl;

        Kernel& kernel = runTime.initialize_kernel();
        compute_kernel(logger, 1.5*dx, particle_neighbours, kernel);

        sid = runTime.get_searchCubeIds();

        compute_density(logger, particle_neighbours, kernel, rho);
        compute_pressure(logger, particle_neighbours, rho,  p);
        compute_pressure_gradient(logger, particle_neighbours, kernel, rho, p, dp);
        compute_dnu(logger, particle_neighbours, kernel, u, rho, dnu, dx);
        compute_du(logger, particle_neighbours, kernel, du, u, rho, dnu, dp, dx);
        compute_u(logger, type, du, u, dt);
        update_pos(logger, points, u, dt);
        // if (runTime.get_timeStep() % 10 == 0) {
            runTime.write_to_disk();
        // }
        particle_neighbours = runTime.update_particle_neighbours(points);
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

