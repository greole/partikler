#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <assert.h>


// #include <future>
// #include <execution>
#include <algorithm>
#include <math.h>    // std::min
// #include <parallel/settings.h>
#include <omp.h>

#include "SPHDatastructures.hpp"
#include "CGALHelper.hpp"
#include "SPHio.hpp"
#include "SPHCore.hpp"
#include "SPHModels.hpp"
#include "computes.cpp"

#include "include/particle_helper.hpp"

/* file IO */
#include <sys/stat.h>
#include <pthread.h>
#include <math.h>
#include <x86intrin.h>

// getopts
#include <getopt.h>


void print_help() {
  std::cout <<
    "--dx <val>:           Average particle distance\n"
    "--os <val>:           Oversubscribe particles by val %\n"
    "--dt <val>:           Max timestep\n"
    "--ts <val>:           Number of timesteps\n"
    "--nl <val>:           Number of extrusion layer\n"
    "--nw <val>:           Writeout frequency\n"
    "--stl <val>:          STL filename\n"
    "--help:              Show help\n";
    exit(1);
}

struct SimArgs {
  float dx;
  float dt;
  float os;
  int ts;
  int nl;
  int nw;
  std::string stl;
};

void print_sim_args(SimArgs s) {
  std::cout
    << "dx: " << s.dx
    << "\nos: " << s.os
    << "\ndt: " << s.dt
    << "\nts: " << s.ts
    << "\nnl: " << s.nl
    << "\nnw: " << s.nw
    << "\nstl: " << s.stl
    << std::endl;
}

SimArgs process_args (int argc, char** argv)
{
    const char* const short_opts = "x:o:t:n:l:f:s:h";
    const option long_opts[] = {
            {"dx", required_argument, nullptr, 'x'},
            {"os", required_argument, nullptr, 'o'},
            {"dt", required_argument, nullptr, 't'},
            {"ts", required_argument, nullptr, 'n'},
            {"nl", required_argument, nullptr, 'l'},
            {"nw", required_argument, nullptr, 'f'},
            {"stl", required_argument, nullptr, 's'},
            {nullptr, no_argument, nullptr, 0}
    };

    // Default args
    float dx = 1.0;
    float dt = 1.0;
    float os = 0.0;
    int ts = 100;
    int nl = 1;
    int nw = 1;
    std::string stl = "in.stl";

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'x':
            // num = std::stoi(optarg);
            dx = std::stof(optarg);
            break;

        case 'o':
          os = std::stof(optarg);
          break;


        case 't':
            dt = std::stof(optarg);
            break;

        case 'n':
            ts = std::stoi(optarg);
            break;

        case 'l':
            nl = std::stoi(optarg);
            break;

        case 'f':
          nw = std::stoi(optarg);
          break;

        case 's':
          stl = std::string(optarg);
          break;

        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            print_help();
            break;
        }
    }
    return {dx, dt, os, ts, nl, nw, stl};
}


// TODO move to model part
SPHPointField generate_boundary_particles(
        SimArgs arg
){
    print_sim_args(arg);
    Logger logger {1};
    RunTime runTime {logger, false};
    float dx = arg.dx;

    // TODO move this to a separate lib
    logger.info_begin() << "Reading input file: " <<  arg.stl;

    std::ifstream istream(arg.stl);

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
    float kernel_relaxation = 1.0;
    float noise_relaxation = 1.0;

    runTime.set_dict("dx", kernel_relaxation*dx);
    runTime.set_dict("n_timesteps", arg.ts);
    runTime.set_dict("write_freq", arg.nw);

    const size_t number_of_facets = std::distance(
            polyhedron.facets_begin(),
            polyhedron.facets_end());

    logger.info() << "Number of facets: " << number_of_facets;

    std::vector<Point>   points;
    std::vector<Facet_handle> initial_facets{};
    initial_facets.reserve(points.size());
    // std::vector<CGALVector>  facet_vectors(number_of_facets);

    Generate_Points_at_Facets gpf(dx, arg.os, points, initial_facets);

    std::for_each(
        polyhedron.facets_begin(),
        polyhedron.facets_end(),
        gpf
    );

    // SPHModelGraph modelGraph("Main");

    // SPHModelGraph setup("Setup");


    SPHModelFactory::createInstance("SPHModelGraph");

    // modelGraph.push_back(setup);

    // modelGraph.execute();

    runTime.set_facets(initial_facets);

    std::vector<Facet_handle>& facets = runTime.get_facets();

    const size_t n_points = points.size();

    logger.info() << "Number of particles: " <<  n_points;

    Vector zeroVec = {0, 0, 0};

    SortedNeighbours &particle_neighbours =
      runTime.initialize_particle_neighbours(points);

    auto sIdxs = runTime.get_sorting_idxs();
    // reorder_vector(sIdxs, facets);

    SPHPointField &pos = runTime.get_particle_positions();

    Kernel &kernel = runTime.initialize_kernel();

    SPHFloatField &rho = runTime.create_uniform_field(0.0, "rho");
    SPHFloatField &p = runTime.create_uniform_field(1.0e5,  "p");

    SPHFloatField &nu = runTime.create_uniform_field(1e4,  "mu");

    SPHVectorField &dp =
      runTime.create_uniform_field(zeroVec, "dp", {"dpx", "dpy", "dpz"});

    SPHVectorField &u =
      runTime.create_uniform_field(zeroVec, "Vel", {"U", "V", "W"});

    SPHVectorField &du =
      runTime.create_uniform_field(zeroVec, "dVel", {"dU", "dV", "dW"});

    SPHVectorField &dnu =
      runTime.create_uniform_field(zeroVec, "dnu", {"dnux", "dnuy", "dnuz"});

    SPHVectorField &f =
      runTime.create_uniform_field(zeroVec, "f", {"fx", "fy", "fz"});

    // surface slide particle have type 2
    SPHIntField &type = runTime.create_uniform_int_field(2, "type");

    SPHSizeTField &idx = runTime.create_idx_field();

    idx.reorder_field(runTime.get_sorting_idxs());
    type.reorder_field(runTime.get_sorting_idxs());

    compute_kernel(
        logger,
        dx * kernel_relaxation,
        runTime.get_particle_positions(),
        particle_neighbours,
        facets,
        kernel);

    float dt;
    dt = 1e-32;
    while (!runTime.end()) {

        runTime.print_timestep();

        // dt = min(dtmax, (float)(1.0e-12*pow(1.05, (float) runTime.get_timestep())));
        // dt = min(arg.dt, (float)(arg.dt/1000.0 + arg.dt/100000.0*(float)runTime.get_timestep()));
        std::cout << "CURRENT TIMESTEP " << dt << std::endl;

        compute_density(logger, particle_neighbours, kernel, rho);

        compute_pressure(logger, particle_neighbours, rho, p);

        compute_pressure_gradient(
            logger, particle_neighbours, kernel, rho, p, dp);

        compute_dnu(
            logger,
            particle_neighbours,
            kernel,
            u,
            rho,
            nu,
            runTime.get_particle_positions(),
            dnu);

        compute_forces(
            logger,
            particle_neighbours,
            runTime.get_particle_positions(),
            facets,
            dx*kernel_relaxation,
            f);

        compute_du(logger, particle_neighbours, kernel, u, f, rho, p, dnu, dp, du);

        // limit_dt( du, dx, dt);

        compute_u(logger, du, dt, u);

        noise_relaxation *= 0.999;
        add_random_noise(logger, u, dx*0.01*noise_relaxation);

        dt = update_pos(logger, u, dt, dx, points, type, facets, idx, pos);

        runTime.write_to_disk();

        // Nullify u
        // u.set(zeroVec);
        u = u/1.01;

        runTime.update_neighbours();

        // kernel_relaxation = max(1.0, kernel_relaxation*0.95);

        runTime.set_dict("dx", dx*kernel_relaxation);

        // du.reorder_field(runTime.get_sorting_idxs());
        u.reorder_field(runTime.get_sorting_idxs());
        idx.reorder_field(runTime.get_sorting_idxs());
        type.reorder_field(runTime.get_sorting_idxs());

        logger.info_begin() << "reordering particle info";

        auto sIdxs = runTime.get_sorting_idxs();
        logger.info_end();

        compute_kernel(
                       logger,
                       dx * kernel_relaxation,
                       runTime.get_particle_positions(),
                       particle_neighbours,
                       facets,
                       kernel);

        runTime++;
    }

    return pos;
}

int main(int argc, char* argv[]) {
    // Step 1 generate boundary particles
    auto boundary_particles = generate_boundary_particles(process_args(argc, argv));
    // Step 2 fix boundary particles extrude and rerun with Lennard-Jones Potential
}
