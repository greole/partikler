/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

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

// external header
#include "yaml-cpp/yaml.h"
#include <getopt.h>

#include "Object.hpp"
#include "Models.hpp"

// #include "SPHDatastructures.hpp"
// #include "CGALHelper.hpp"

// #include "SPHio.hpp"
// #include "SPHCore.hpp"
// #include "ParticleNeighbours.hpp"

// #include "computes.cpp"

// #include "include/particle_helper.hpp"

// /* file IO */
// #include <sys/stat.h>
// #include <pthread.h>
// #include <math.h>
// #include <x86intrin.h>



void print_help() {
  std::cout <<
    "--config <path/to/SPH.yaml>:             Location of config file\n";
    "--help:                                  Show help\n";
    exit(1);
}

YAML::Node process_args (int argc, char** argv)
{
    const char* const short_opts = "x:o:t:n:l:f:s:h";
    const option long_opts[] = {
            {"config", required_argument, nullptr, 'c'},
            {nullptr, no_argument, nullptr, 0}
    };

    // Default args
    std::string conf;

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {

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

// // TODO move to model part
// PointField generate_boundary_particles(
//     YAML::Node parameter
// ){
//     Logger logger {1};

//     ObjectRegistry obj_reg {};

//     // RunTime runTime {logger, false};

//     int ts = parameter["timesteps"].as<int>();;
//     int nw = parameter["writeout"].as<int>();;

//     obj_reg.create_generic<Generic<TimeInfo>>(
//         "TimeInfo", TimeInfo {1e-32, ts, 1e24});

//     TimeGraph loop {"TimeGraph", YAML::Node(), obj_reg};

//     // Register Models

//     for (auto el: parameter["ModelGraph"]["pre"]) {

//         auto model_namespace = el.first.as<std::string>();
//         auto model_name = el.second["model"].as<std::string>();
//         auto params = el.second;

//         auto model = ModelFactory::createInstance(
//             model_namespace,
//             model_name,
//             model_name,
//             el.second,
//             obj_reg
//             // runTime
//             );

//         loop.push_back_pre(model);
//     }

//     loop.execute_pre();

//     // Part 1:
//     // Distribute points randomly over cell facets/triangles
//     // Returns a vector of initial point packets per facet
//     // float kernel_relaxation = 1.0;
//     float noise_relaxation = 1.0;


//     // surface slide particle have type 2
//     IntField &type = obj_reg.create_field<IntField>("type", 2);
//     SizeTField &idx = obj_reg.create_idx_field();

//     // // Register Models
//     for (auto el: parameter["ModelGraph"]["main"]) {

//         // auto model_name = el.first.as<std::string>();
//         // auto params = el.second;
//         auto model_namespace = el.first.as<std::string>();
//         auto model_name = el.second["model"].as<std::string>();

//         auto model = ModelFactory::createInstance(
//             model_namespace,
//             model_name,
//             model_name,
//             el.second,
//             obj_reg);

//         loop.push_back_main(model);
//     }

//     loop.execute_main();

//     return obj_reg.get_particle_positions();
// }

int main(int argc, char* argv[]) {
    // Step 1 generate boundary particles
    YAML::Node config = process_args(argc, argv);

    ObjectRegistry obj_reg {};

    TimeGraph &timeGraph =
        obj_reg.register_object<TimeGraph>(std::make_unique<TimeGraph>(
            "TimeGraph", config["PROJECT"], obj_reg));

    // main model loop
    for (auto el: config["PROJECT"]["PRE"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();
        auto params = it->second;

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            params,
            obj_reg
            );

        timeGraph.push_back_pre(model);
    }

    timeGraph.execute_pre();

    std::cout << "main" << std::endl;

    for (auto el: config["PROJECT"]["MAIN"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();
        auto params = it->second;

        std::cout << model_namespace << std::endl;
        std::cout << model_name << std::endl;

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            params,
            obj_reg
            );

        timeGraph.push_back_main(model);
    }

    timeGraph.execute_main();
}
