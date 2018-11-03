#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>

// #include <future>
// #include <execution>
#include <algorithm>
// #include <parallel/settings.h>
#include <omp.h>

#include "SPHio.h"
#include "SPHCore.h"
#include "include/particle_helper.h"

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

    RunTime runTime {1};

    // TODO move this to a separate lib
    runTime.info_begin() << "Reading input file: " <<  argv[1];

    std::ifstream istream(argv[1]);

    // read_STL(istream, points, facets, false);
    Polyhedron_builder_from_STL<HalfedgeDS> builder(istream);
    runTime.info_end();

    runTime.info_begin() << "Constructing polyhedron";
    // Create input polyhedron
    Polyhedron polyhedron;
    polyhedron.delegate(builder);
    runTime.info_end();

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    const float dx = atof(argv[3]);

    const size_t number_of_facets = std::distance(
            polyhedron.facets_begin(),
            polyhedron.facets_end());

    runTime.info() << "Number of facets: " << number_of_facets;

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

    runTime.info() << "Number of particles: " <<  n_points;

    // TODO reimplement layer extrusion

    // SPH steps
    // compute distance to nearest neighbour points
    // generates 2d array distances [npoints][default_neigbours]

    // initialise search cubes
    // an array of particle ids and cgal circulator like iterator
    // search_cube[id][particle_id]
    // supports search_cube.cube(id).begin();
    // search_cubes
    //
    SearchCubeTree search_cubes (points, 3*dx, RunTime(1));

    const ParticleNeighbours& particle_neighbours =
        search_cubes.get_particle_distances();

    // compute and store the kernel

    runTime.info_begin() << "Initialising Kernel";
    Kernel kernel;
    kernel.W    = std::move(std::vector<float>  (particle_neighbours.origId.size(), 0));
    kernel.dWdx = std::move(std::vector<Vector> (particle_neighbours.origId.size(), {0, 0, 0}));
    runTime.info_end();

    compute_kernel(runTime, 3*dx, particle_neighbours, kernel);

    runTime.info_begin() << "Initialising particle densities";
    std::vector<float> densities (n_points, 0.0);
    runTime.info_end();



    runTime.info_begin() << "Initialising particle velocities";
    std::vector<float>  nu (particle_neighbours.origId.size(), 0);
    std::vector<Vector> velocities (n_points, {0, 0, 0});
    std::vector<Vector> du (n_points, {0, 0, 0});
    runTime.info_end();

    for (int i=0; i<20; i++) {
        compute_density(runTime, particle_neighbours, kernel, densities);
        compute_nu(runTime, particle_neighbours, velocities, densities, nu, dx);
        compute_du(runTime, particle_neighbours, kernel, du, velocities, densities, nu, dx);
        compute_u(runTime, particle_neighbours, kernel, du, velocities, densities, nu, dx, 1000000);
        writeData_SPH("daten", i, points, densities, nu, du, velocities);
    }
    return 0;
}
