#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>

// #include <future>
// #include <execution>
#include <algorithm>
// #include <parallel/settings.h>
#include <omp.h>

#include "SPHDatastructures.hpp"
#include "SPHio.hpp"
#include "SPHCore.hpp"
#include "computes.cpp"

#include "include/particle_helper.hpp"

/* file IO */
#include <sys/stat.h>
#include <pthread.h>
#include <math.h>
#include <x86intrin.h>



void generate_boundary_particles(
        bool verbose,
        bool write_substeps,
        std::string stl_filename,
        float dx,
        float dt,
        size_t nts,
        int n_extrusion_layer
){
    Logger logger {1};
    RunTime runTime {logger, false};

    // TODO move this to a separate lib
    logger.info_begin() << "Reading input file: " <<  stl_filename;

    std::ifstream istream(stl_filename);

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
    runTime.set_dict("dx", dx);
    runTime.set_dict("n_timesteps", nts);

    const size_t number_of_facets = std::distance(
            polyhedron.facets_begin(),
            polyhedron.facets_end());

    logger.info() << "Number of facets: " << number_of_facets;

    //TODO std::vector:reserve
    std::vector<Point>   points;
    std::vector<size_t>  number_points_facet(number_of_facets);
    std::vector<Facet_handle>   initial_facets(number_of_facets);
    std::vector<CGALVector>  facet_vectors(number_of_facets);

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

    // Generate "Free" Particles by cloning original particles 
    // and shifting slightly outwards each layer has max distances
    // to original particle

    // Make a copy of original sorted particles to restore

    // SPHPointField positions = runTime.set_particle_positions(points);

    FixedDistanceParticles extruded_points =
        create_extruded_points(points, n_extrusion_layer, 0.5*dx, gpf);

    // Generate_Points_at_Edges  pae(dx, extruded_points);

    // std::for_each(
    //               polyhedron.facets_begin(),
    //               polyhedron.facets_end(),
    //               pae
    //               );

    Vector zeroVec = {0, 0, 0};

    SortedNeighbours &particle_neighbours =
      runTime.initialize_particle_neighbours(extruded_points.points);

    SPHPointField &pos = runTime.get_particle_positions();

    Kernel &kernel = runTime.initialize_kernel();

    SPHFloatField &rho = runTime.create_uniform_field(0.0, "rho");
    SPHFloatField &p = runTime.create_uniform_field(1.0e5,  "p");
    SPHVectorField &dp =
      runTime.create_uniform_field(zeroVec, "dp", {"dpx", "dpy", "dpz"});
    SPHVectorField &u =
      runTime.create_uniform_field(zeroVec, "Vel", {"U", "V", "W"});

    SPHVectorField &du =
      runTime.create_uniform_field(zeroVec, "dVel", {"dU", "dV", "dW"});

    SPHSizeTField &idx = runTime.create_idx_field();

    compute_kernel(
        logger,
        2.8 * dx,
        runTime.get_particle_positions(),
        particle_neighbours,
        kernel);

    // compute "rubber band" forces
    // regular sph step but distance to original particle is limited
    auto sIdxs = runTime.get_sorting_idxs();

    // TODO reorder(sIdxs, extruded_points) via olverload
    reorder_vector(sIdxs, extruded_points.fixId);//
    reorder_vector(sIdxs, extruded_points.maxDx);//
    reorder_vector(sIdxs, extruded_points.dir);//
    reorder_vector(sIdxs, extruded_points.mType);//
    reorder_vector(sIdxs, extruded_points.facets);//

    // timestep loop
    while (!runTime.end()) {

        logger.info() << "Timestep: ";
        std::cout << runTime.get_timestep() << std::endl;

        compute_density(logger, particle_neighbours, kernel, rho);

        compute_pressure(logger, particle_neighbours, rho, p);

        compute_pressure_gradient(
            logger, particle_neighbours, kernel, rho, p, dp);

        compute_du(logger, particle_neighbours, kernel, u, rho, p, dp, du);

        compute_u(logger, du, dt, u);

        update_pos(logger, u, dt, points, extruded_points, initial_facets, idx, pos);

        runTime.write_to_disk();

        runTime.update_neighbours();

        compute_kernel(
                       logger,
                       2.8 * dx,
                       runTime.get_particle_positions(),
                       particle_neighbours,
                       kernel);

        // du.reorder_field(runTime.get_sorting_idxs());
        u.reorder_field(runTime.get_sorting_idxs());
        idx.reorder_field(runTime.get_sorting_idxs());

        logger.info_begin() << "reordering particle info";
        auto sIdxs = runTime.get_sorting_idxs();

        reorder_vector(sIdxs, extruded_points.fixId);//
        reorder_vector(sIdxs, extruded_points.maxDx);//
        reorder_vector(sIdxs, extruded_points.dir);//
        reorder_vector(sIdxs, extruded_points.mType);//
        reorder_vector(sIdxs, extruded_points.facets);//
        logger.info_end();

        runTime++;
    }
}

int main(int argc, char* argv[]) {

    const float dx = atof(argv[2]);
    const float dt = atof(argv[3]);
    const size_t n_timesteps = atof(argv[4]);
    const int layer = atof(argv[5]);
    generate_boundary_particles(true, true, argv[1], dx, dt, n_timesteps, layer);
}
