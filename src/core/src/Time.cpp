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
#include "Time.hpp"

TimeGraph::TimeGraph(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      init_(ModelGraph("pre", parameter, objReg)),
      main_(ModelGraph("main", parameter, objReg)),
      post_(ModelGraph("post", parameter, objReg)), current_timestep_(0),
      current_time_(0), name_(read_coeff<std::string>("name")),
      endTime_(read_or_default_coeff<Scalar>("endTime", -1.0)),
      deltaT_(read_or_default_coeff<Scalar>("deltaT", -1.0)),
      max_deltaT_(read_or_default_coeff<Scalar>("max_deltaT", -1.0)),
      iterations_(read_or_default_coeff<int>("iterations", 0)) {
    if (deltaT_ < 0) {
        iter_mode = true;
    }
}

void TimeGraph::execute() {
    execute_pre();
    execute_main();
    execute_post();
}

REGISTER_DEF_TYPE(CORE, TimeGraph);
