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

#include "Gravity.hpp"

#include "Time.hpp"

Gravity::Gravity(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : VectorFieldEquation(
          model_name,
          parameter,
          objReg,
          objReg.create_field<VectorField>(
              "g", zero<VectorField::value_type>::val, {"gx", "gy", "gz"})),
      gravity_(read_or_default_vec3("g", {0, 0, 0})) {

    for (auto &el : f_) {
        el = gravity_;
    }
}

void Gravity::execute() {

    log().info_begin() << "Computing gravity";

    iteration_ = time_.get_current_timestep();

    log().info_end();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Gravity);
