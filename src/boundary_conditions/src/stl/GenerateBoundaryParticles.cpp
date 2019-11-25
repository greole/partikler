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



#include "GenerateBoundaryParticles.hpp"


GenerateBoundaryParticles::GenerateBoundaryParticles(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : Model(model_name, parameter, objReg),
      boundaryIds_(objReg.create_field<IntField>("boundary")),
      typeIds_(objReg.create_field<IntField>("type")),
      idx_(objReg.create_field<SizeTField>("idx")),
      pos_(objReg.create_field<PointField>("Pos")),
      ts_(read_coeff<int>("iterations")),
      nw_(read_or_default_coeff<int>("writeout", -1))
{};

void GenerateBoundaryParticles::execute() {

    Logger logger {1};

    get_objReg().create_generic<Generic<TimeInfo>>(
        "TimeInfo", TimeInfo {1e-32, ts_, 1e24});

    TimeGraph loop {"TimeGraph", YAML::Node(), get_objReg()};

    // Register Models
    for (auto el: get_parameter()) {
        std::cout << el << std::endl;
    }

    for (auto el: get_parameter()["ModelGraph"]["pre"]) {

        auto model_namespace = el.first.as<std::string>();
        auto model_name = el.second["model"].as<std::string>();
        auto params = el.second;

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            el.second,
            get_objReg()
            );

        loop.push_back_pre(model);
    }

    loop.execute_pre();

    std::cout << pos_.size() << std::endl;

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    // float kernel_relaxation = 1.0;
    float noise_relaxation = 1.0;


    // surface slide particle have type 2
    // IntField &type = get_objReg().create_field<IntField>("type", 2);
    // SizeTField &idx = obj_reg.create_idx_field();

    // // Register Models
    for (auto el: get_parameter()["ModelGraph"]["main"]) {

        // auto model_name = el.first.as<std::string>();
        // auto params = el.second;
        auto model_namespace = el.first.as<std::string>();
        auto model_name = el.second["model"].as<std::string>();

        auto model = ModelFactory::createInstance(
            model_namespace,
            model_name,
            model_name,
            el.second,
            get_objReg());

        loop.push_back_main(model);
    }

    loop.execute_main();

    // update idx and typeIds
};

REGISTER_DEF_TYPE(BOUNDARY, GenerateBoundaryParticles);
