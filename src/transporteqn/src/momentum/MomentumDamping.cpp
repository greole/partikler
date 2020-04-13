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

#include "MomentumDamping.hpp"
#include "Equation.hpp"

#include "Time.hpp"

MomentumDamping::MomentumDamping(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : VectorFieldEquation(
          "MomentumDamping", parameter, objReg, objReg.velocity()),
      nu_(read_or_default_coeff<Scalar>("nu", 1.0)) {}

void MomentumDamping::execute() {

    log().info_begin() << "Computing du/dt";

    for (auto &el : f_) {
        el = nu_ * el;
    }

    iteration_ = time_.get_current_timestep();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, MomentumDamping);
