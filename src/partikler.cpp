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
#include "ParticleNeighbours.hpp"

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
    // "--dx <val>:           Average particle distance\n"
    // "--os <val>:           Oversubscribe particles by val %\n"
    // "--dt <val>:           Max timestep\n"
    // "--ts <val>:           Number of timesteps\n"
    // "--nl <val>:           Number of extrusion layer\n"
    // "--nw <val>:           Writeout frequency\n"
    // "--stl <val>:          STL filename\n"
    "--config <path/to/SPH.yaml>:             Location of config file\n";
    "--help:                                  Show help\n";
    exit(1);
}

// struct SimArgs {
//   float dx;
//   float dt;
//   float os;
//   int ts;
//   int nl;
//   int nw;
//   std::string stl;
// };

// void print_sim_args(SimArgs s) {
//   std::cout
//     << "dx: " << s.dx
//     << "\nos: " << s.os
//     << "\ndt: " << s.dt
//     << "\nts: " << s.ts
//     << "\nnl: " << s.nl
//     << "\nnw: " << s.nw
//     << "\nstl: " << s.stl
//     << std::endl;
// }

YAML::Node process_args (int argc, char** argv)
{
    const char* const short_opts = "x:o:t:n:l:f:s:h";
    const option long_opts[] = {
            {"config", required_argument, nullptr, 'c'},
            // {"os", required_argument, nullptr, 'o'},
            // {"dt", required_argument, nullptr, 't'},
            // {"ts", required_argument, nullptr, 'n'},
            // {"nl", required_argument, nullptr, 'l'},
            // {"nw", required_argument, nullptr, 'f'},
            // {"stl", required_argument, nullptr, 's'},
            {nullptr, no_argument, nullptr, 0}
    };

    // Default args
    // float dx = 1.0;
    // float dt = 1.0;
    // float os = 0.0;
    // int ts = 100;
    // int nl = 1;
    // int nw = 1;
    std::string conf;

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        // case 'x':
        //     // num = std::stoi(optarg);
        //     dx = std::stof(optarg);
        //     break;

        // case 'o':
        //   os = std::stof(optarg);
        //   break;


        // case 't':
        //     dt = std::stof(optarg);
        //     break;

        // case 'n':
        //     ts = std::stoi(optarg);
        //     break;

        // case 'l':
        //     nl = std::stoi(optarg);
        //     break;

        // case 'f':
        //   nw = std::stoi(optarg);
        //   break;

        case 'c':
          conf = std::string(optarg);
          break;

        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            print_help();
            break;
        }
    }

    return YAML::LoadFile(conf);
}


// TODO move to model part
SPHPointField generate_boundary_particles(
    YAML::Node parameter
){
    Logger logger {1};

    SPHObjectRegistry obj_reg {};

    RunTime runTime {logger, false};

    float dx = parameter["ParticleGenerator"]["dx"].as<float>();
    int ts = parameter["timesteps"].as<int>();;
    int nw = parameter["writeout"].as<int>();;

    runTime.create_field<SPHPointField>("Pos");

    runTime.create_field<SPHField<Facet_handle>>("facets");


    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    // float kernel_relaxation = 1.0;
    float noise_relaxation = 1.0;

    // runTime.set_dict("dx", kernel_relaxation*dx);
    runTime.set_dict("dx", dx);
    runTime.set_dict("n_timesteps", ts);
    runTime.set_dict("write_freq",  nw);

    runTime.create_generic<SPHGeneric<TimeInfo>>(
        "TimeInfo", TimeInfo {1e-32, ts, 1e24});

    runTime.create_generic<SPHGeneric<SearchCubeDomain>>("search_cube_domain");

    auto main = SPHModelFactory::createInstance(
        "CORE",
        "SPHModelGraph",
        "main",
        YAML::Node(),
        runTime);

    auto pre = SPHModelFactory::createInstance(
        "CORE",
        "SPHModelGraph",
        "pre",
        YAML::Node(),
        runTime);

    auto reader = SPHModelFactory::createInstance(
        "READER",
        "SPHSTLReader",
        "SPHSTLReader",
        parameter["STLReader"],
        runTime);

    auto generator = SPHModelFactory::createInstance(
        "GENERATOR",
        "SPHParticleGenerator",
        "ParticleGenerator",
        parameter["ParticleGenerator"],
        runTime);

    main->execute();
    reader->execute();
    generator->execute();

    auto kernel_field = runTime.create_field<SPHFloatField>( "KernelW", 0.0);

    kernel_field.set_reorder(false);

    runTime.create_field<SPHField<VectorPair>>("KerneldWdx");
    runTime.create_field<SPHField<NeighbourPair>>("neighbour_pairs");
    runTime.create_field<SPHField<STLSurfaceDist>>("surface_dist");
    runTime.create_field<SPHField<SearchCube>>("search_cubes");
    runTime.create_field<SPHSizeTField>("sorting_idxs");

    auto kernel = SPHModelFactory::createInstance(
        "KERNEL",
        "STLWendland2D",
        "Kernel",
        parameter["STLWendland2D"],
        runTime);


    auto neighbours = SPHModelFactory::createInstance(
        "PARTICLENEIGHBOURS",
        "SPHSTLParticleNeighbours",
        "Neighbours",
        parameter["SPHSTLParticleNeighbours"],
        runTime);

    neighbours->execute();

    kernel->execute();

    runTime.update();

    Vector zeroVec = {0, 0, 0};

    runTime.create_field<SPHFloatField>("rho", 0.0);
    runTime.create_field<SPHFloatField>("p", 1.0e5);
    runTime.create_field<SPHFloatField>("mu", 1e4);
    runTime.create_field<SPHVectorField>("dp", zeroVec, {"dpx", "dpy", "dpz"});
    runTime.create_field<SPHVectorField>("du", zeroVec, {"dU", "dV", "dW"});
    runTime.create_field<SPHVectorField>("dnu", zeroVec,  {"dnux", "dnuy", "dnuz"});
    runTime.create_field<SPHVectorField>("f", zeroVec, {"fx", "fy", "fz"});
    auto & u = runTime.create_field<SPHVectorField>("u", zeroVec, {"U", "V", "W"});

    // surface slide particle have type 2
    SPHIntField &type = runTime.create_field<SPHIntField>("type", 2);

    SPHSizeTField &idx = runTime.create_idx_field();

    auto conti = SPHModelFactory::createInstance(
        "TRANSPORTEQN",
        "Conti",
        "Conti",
        parameter["Conti"],
        runTime);


    auto pressure = SPHModelFactory::createInstance(
        "TRANSPORTEQN",
        "Pressure",
        "Pressure",
        parameter["Pressure"],
        runTime);

    auto momentum = SPHModelFactory::createInstance(
        "TRANSPORTEQN",
        "Momentum",
        "Momentum",
        parameter["Momentum"],
        runTime);

    auto visc = SPHModelFactory::createInstance(
        "TRANSPORTEQN",
        "Viscosity",
        "Viscosity",
        parameter["Viscosity"],
        runTime);

    auto pos = SPHModelFactory::createInstance(
        "TRANSPORTEQN",
        "STLPosIntegrator",
        "STLPosIntegrator",
        parameter["STLPosIntegrator"],
        runTime);

    kernel_field.set_reorder(false);

    conti->execute();

    pressure->execute();

    visc->execute();

    momentum->execute();

    pos->execute();

    runTime.write_to_disk();

    float dt;
    dt = 1e-32;

    runTime++;
    while (!runTime.end()) {

        runTime.print_timestep();

        neighbours->execute();

        kernel->execute();

        conti->execute();

        pressure->execute();

        visc->execute();

        momentum->execute();

        noise_relaxation *= 0.999;
        add_random_noise(logger, u, dx*0.01*noise_relaxation);

        pos->execute();

        runTime.write_to_disk();
        runTime++;
    }

    return runTime.get_particle_positions();
}

int main(int argc, char* argv[]) {
    // Step 1 generate boundary particles
    YAML::Node config = process_args(argc, argv);
    auto boundary_particles =
        generate_boundary_particles(config["GenerateBoundaryParticles"]);
}
