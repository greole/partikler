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

#include <execinfo.h>
#include <getopt.h> // for getopt_long, no_argument
#include <iostream> // for endl, operator<<, cout
#include <memory>   // for make_unique, allocator
#include <signal.h>
#include <stdlib.h>
#include <string> // for string, operator<<

#include "Models.hpp"         // for ModelFactory
#include "Object.hpp"         // for GenericType
#include "ObjectRegistry.hpp" // for FieldIdMap, ObjectReg...
#include "Time.hpp"
#include "yaml-cpp/node/detail/iterator.h"     // for iterator_base, iterat...
#include "yaml-cpp/node/detail/iterator_fwd.h" // for const_iterator
#include "yaml-cpp/node/impl.h"                // for Node::operator[], Nod...
#include "yaml-cpp/node/iterator.h"            // for iterator_value
#include "yaml-cpp/node/node.h"                // for Node
#include "yaml-cpp/node/parse.h"               // for LoadFile

#include "Scalar.hpp"

// // Currently under ubuntu the gcc compiler optimises
// // the register to call to the model register away
// // if it doesn't see a call anywere in the executable.
// // Hence, for no the register calls are made manually
// // until a better solution is found
// #ifdef WITH_GNU
// #include "GenerateBoundaryParticles.hpp"

// REGISTER_DEF_TYPE(BOUNDARY, GenerateBoundaryParticles);
// #include "Wendland2D.hpp"

// REGISTER_DEF_TYPE(KERNEL, STLWendland2D);
// #include "ParticleNeighbours.hpp"
// #include "STLParticleNeighbours.hpp"

// REGISTER_DEF_TYPE(PARTICLENEIGHBOURS, SPHSTLParticleNeighbours);
// REGISTER_DEF_TYPE(PARTICLENEIGHBOURS, SPHParticleNeighbours);
// #include "SortParticles.hpp"

// REGISTER_DEF_TYPE(SORTING, CountingSortParticles);
// #include "Conti.hpp"

// REGISTER_DEF_TYPE(TRANSPORTEQN, Conti);
// #include "Momentum.hpp"

// REGISTER_DEF_TYPE(TRANSPORTEQN, Momentum);
// #include "Bonet.hpp"
// #include "Cole.hpp"
// #include "Solenthaler.hpp"

// REGISTER_DEF_TYPE(TRANSPORTEQN, Cole);
// #include "Szewc.hpp"

// REGISTER_DEF_TYPE(VISCOSITY, Szewc);
// #include "CGALParticleGenerator.hpp"

// REGISTER_DEF_TYPE(READER, SPHSTLReader);
// REGISTER_DEF_TYPE(GENERATOR, SPHParticleGenerator);
// #include "STLPosIntegrator.hpp"

// REGISTER_DEF_TYPE(TRANSPORTEQN, STLPosIntegrator);
// #include "SuperSPHWriter.hpp"

// REGISTER_DEF_TYPE(EXPORT, SuperSPHWriter);

// #include "HDF5Writer.hpp"
// REGISTER_DEF_TYPE(EXPORT, HDF5Writer);
// #include "CreateFields.hpp"

// REGISTER_DEF_TYPE(FIELDS, InitFields);
// #include "FixedValue.hpp"

// REGISTER_DEF_TYPE(BOUNDARY, FixedValue);

// #include "Cubiod.hpp"
// REGISTER_DEF_TYPE(FIELDS, InitShape);

// #endif

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

void print_help() {
    std::cout
        << "--config <path/to/SPH.yaml>:             Location of config file\n";
    "--help:                                  Show help\n";
    exit(1);
}

YAML::Node process_args(int argc, char **argv) {
    const char *const short_opts = "x:o:t:n:l:f:s:h";
    const option long_opts[] = {{"config", required_argument, nullptr, 'c'},
                                {nullptr, no_argument, nullptr, 0}};

    // Default args
    std::string conf;

    while (true) {
        const auto opt =
            getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt) break;

        switch (opt) {

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

int main(int argc, char *argv[]) {
    signal(SIGSEGV, handler); // install our handler
    // Step 1 generate boundary particles
    YAML::Node config = process_args(argc, argv);

    ObjectRegistry obj_reg {};

    // Register Materials

    MaterialMap &mm = obj_reg.register_object<MaterialMap>(
        std::make_unique<MaterialMap>("MaterialMap", GenericType));

    for (auto el : config["PROJECT"]["MATERIALS"]) {
        YAML::const_iterator it = el.begin();
        std::string material_name = it->first.as<std::string>();
        std::string material_type = it->second["type"].as<std::string>();
        Scalar mp = it->second["mp"].as<Scalar>();
        Scalar rho = it->second["rho"].as<Scalar>();

        Material m {material_name, material_type, mp, rho};

        mm.insert(material_name, m);
    }

    TimeGraph &timeGraph = obj_reg.register_object<TimeGraph>(
        std::make_unique<TimeGraph>("TimeGraph", config["PROJECT"], obj_reg));

    // Upon creation every "block" of particles gets its own id
    // so that boundary/solid particles can be identified.
    // Additionally every id maps to a material type.
    FieldIdMap &fieldIdMap = obj_reg.register_object<FieldIdMap>(
        std::make_unique<FieldIdMap>("FieldIdMap", GenericType));

    // main model loop
    for (auto el : config["PROJECT"]["PRE"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();
        auto params = it->second;

        auto model = ModelFactory::createInstance(
            model_namespace, model_name, model_name, params, obj_reg);

        timeGraph.push_back_pre(model);
    }

    timeGraph.execute_pre();

    auto &idx = obj_reg.create_field<ScalarField>("idx", 0.0);
    for (size_t i = 0; i < idx.size(); i++) {
        idx[i] = (Scalar)i;
    }

    for (auto el : config["PROJECT"]["MAIN"]) {

        YAML::const_iterator it = el.begin();
        auto model_namespace = it->first.as<std::string>();
        auto model_name = it->second["model"].as<std::string>();
        auto params = it->second;

        std::cout << model_namespace << std::endl;
        std::cout << model_name << std::endl;

        auto model = ModelFactory::createInstance(
            model_namespace, model_name, model_name, params, obj_reg);

        timeGraph.push_back_main(model);
    }

    timeGraph.execute_main();
}
