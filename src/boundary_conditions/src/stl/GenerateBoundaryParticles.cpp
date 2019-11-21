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
    const std::string &model_name, YAML::Node parameter, RunTime &runTime)

    : SPHModel(model_name, parameter, runTime)
{};

void GenerateBoundaryParticles::execute() {

    int ts = get_params()["timesteps"].as<int>();;
    int nw = get_params()["writeout"].as<int>();;

    TimeGraph loop {"TimeGraph", YAML::Node(), get_runTime()};

    // Register Models

    for (auto el: get_params()["ModelGraph"]["pre"]) {

        auto model_name = el.first.as<std::string>();
        auto params = el.second;

        auto model = SPHModelFactory::createInstance(
            el.second["type"].as<std::string>(),
            model_name,
            model_name,
            el.second,
            get_runTime());

        loop.push_back_pre(model);
    }

    loop.execute_pre();

    // Part 1:
    // Distribute points randomly over cell facets/triangles
    // Returns a vector of initial point packets per facet
    // float kernel_relaxation = 1.0;
    float noise_relaxation = 1.0;

    get_runTime().set_dict("n_timesteps", ts);
    get_runTime().set_dict("write_freq",  nw);

    get_runTime().create_generic<SPHGeneric<TimeInfo>>(
        "TimeInfo", TimeInfo {1e-32, ts, 1e24});

    // surface slide particle have type 2
    SPHIntField &type = get_runTime().create_field<SPHIntField>("type", 2);

    SPHSizeTField &idx = get_runTime().create_idx_field();

    // Register Models
    for (auto el: get_params()["ModelGraph"]["main"]) {

        auto model_name = el.first.as<std::string>();
        auto params = el.second;

        auto model = SPHModelFactory::createInstance(
            el.second["type"].as<std::string>(),
            model_name,
            model_name,
            el.second,
            get_runTime());

        loop.push_back_main(model);
    }

    loop.execute_main();

};

REGISTER_DEF_TYPE(MODEL, GenerateBoundaryParticles);
