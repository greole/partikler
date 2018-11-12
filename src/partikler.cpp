#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>

// #include <future>
// #include <execution>
#include <algorithm>
// #include <parallel/settings.h>
#include <omp.h>

#include "SPHDatastructures.h"
// #include "SPHio.h"
// #include "SPHCore.h"
// #include "include/particle_helper.h"

/* file IO */
#include <sys/stat.h>
#include <pthread.h>
#include <math.h>
#include <x86intrin.h>


// void compute_kernel_sse(
//         RunTime runTime,
//         const float h,
//         const ParticleNeighbours& particle_neighbours,
//         Kernel& kernel
//         )
// {
//     runTime.info_begin() << "Computing Kernel";
//     const size_t size {particle_neighbours.origId.size()};
//     const float invh = 1/h;
//     const float W_fak2  = 21. / (256. * M_PI * h * h * h);
//     const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);
//
//     for (size_t pid = 0; pid < size; pid+=4) {
//
//         // float : 32bit = 4 bytes
//         // load 4 floats into 128 sse register
//         __m128 foo {1.2, 1.0, 2.1, 2.0};
//         //
//         const float len = particle_neighbours.squared_length[pid];
//         const float q {len * invh};
//
//         // if (q > 2.) {
//         //     std::cout << "[DEBUG] outside kernel radius" << std::endl;
//         //     kernel.W[pid] = 0.0;
//         //     kernel.dWdx[pid] = Vector {0.0, 0.0, 0.0};
//         //     continue;
//         // }
//
//         const float q3 = (q - 2.);
//         const float qfac2 = q3*q3;
//         const float qfac4 = qfac2*qfac2;
//
//         float q2 = 2.*q;
//         q2 += 1.;
//
//         kernel.W[pid] = qfac4*q2*W_fak2;
//
//         const float prefact = 10. * qfac2 * q * dW_fak2;
//         kernel.dWdx[pid] = particle_neighbours.normalised_distances[pid]*prefact;
//     }
//
//     runTime.info_end();
// };



int main(int argc, char* argv[]) {

    /*
    Logger logger {1};
    RunTime runTime {logger, false};

    // runTime.setSolverDict();

    // TODO move this to a separate lib
    logger.info_begin() << "Reading input file: " <<  argv[1];

    std::ifstream istream(argv[1]);

    // read_STL(istream, points, facets, false);
    Polyhedron_builder_from_STL<HalfedgeDS> builder(istream);
    logger.info_end();

    logger.info_begin() << "Constructing polyhedron";
    // Create input polyhedron
    Polyhedron polyhedron;
    polyhedron.delegate(builder);
    logger.info_end();

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    const float dx = atof(argv[3]);
    const float dt = atof(argv[4]);
    const size_t n_timesteps = atof(argv[5]);
    runTime.set_dict("dx", dx);
    runTime.set_dict("n_timesteps", n_timesteps);

    const size_t number_of_facets = std::distance(
            polyhedron.facets_begin(),
            polyhedron.facets_end());

    logger.info() << "Number of facets: " << number_of_facets;

    //TODO std::vector:reserve
    std::vector<Point>   points;
    std::vector<size_t>  number_points_facet(number_of_facets);
    std::vector<Facet>   initial_facets(number_of_facets);
    std::vector<Vector>  facet_vectors(number_of_facets);

    Generate_Points_at_Facets gpf(
            dx, points,
            number_points_facet,
            initial_facets,
            facet_vectors
    );

    std::for_each(
        polyhedron.facets_begin(),
        polyhedron.facets_end(),
        gpf
    );

    const size_t n_points = points.size();

    logger.info() << "Number of particles: " <<  n_points;

    SPHPointField positions = runTime.set_particle_positions(points);

    ParticleNeighbours& particle_neighbours = runTime.initialize_particle_neighbours();

    Kernel& kernel = runTime.initialize_kernel();

    compute_kernel(logger, 3*dx, particle_neighbours, kernel);

    std::vector<float>  & rho = runTime.create_uniform_field(1.0      ,"rho").get_field();
    std::vector<float>  &   p = runTime.create_uniform_field(1.01e5   ,  "p").get_field();
    std::vector<Vector> & dp  = runTime.create_uniform_field({0, 0, 0}, "dp").get_field();
    std::vector<Vector> & dnu = runTime.create_uniform_field({0, 0, 0},"dnu").get_field();
    std::vector<Vector> &  du = runTime.create_uniform_field({0, 0, 0}, "du").get_field();
    std::vector<Vector> &   u = runTime.create_uniform_field({0, 0, 0},  "u").get_field();


    for (;!runTime.end(); runTime++) {
        // const float dt = 3e-20;
        runTime.write_to_disk();
        compute_density(logger, particle_neighbours, kernel, rho);
        compute_pressure(logger, particle_neighbours, rho,  p);
        compute_pressure_gradient(logger, particle_neighbours, kernel, rho, p, dp);
        compute_dnu(logger, particle_neighbours, kernel, u, rho, dnu, dx);
        compute_du(logger, particle_neighbours, kernel, du, u, rho, dnu, dp, dx);
        // compute_u(logger, du, u, dt);
        update_pos(logger, points, u, dt);
        particle_neighbours = runTime.update_particle_neighbours(points);
    }
    return 0;
    */
}
