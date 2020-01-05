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

# include "Momentum.hpp"

Momentum::Momentum(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)

    : VectorFieldEquation(model_name, parameter, objReg, objReg.get_object<VectorField>("u")),
      dnu_(objReg.get_object<VectorField>("dnu")),
      dp_(objReg.get_object<VectorField>("dp")),
      u_(objReg.get_object<VectorField>("u")),
      du_(objReg.create_field<VectorField>(
              "du", zero<VectorField::value_type>::val, {"dU", "dV", "dW"})),
      time_(objReg.get_object<TimeGraph>("TimeGraph")) {}

void Momentum::execute() {

    log().info_begin() << "Computing du/dt";

    solve(du_, dnu_ - dp_);

    log().info_end();

    log().info_begin() << "Computing velocity";

    // TODO implement time integrator
     // u_ = u_/1.3;

    // CFL = max(u)*deltaT/dx

    // float maxU = u_.norm().get_max();

    // std::cout << "maxU"
    //           << maxU
    //           << "deltaT" <<  time_().deltaT
    //           << "maxCFL"
    //           << maxU * time_().deltaT /
    //                  std::any_cast<float>(get_runTime().get_dict("dx"))
    //           << std::endl;

    // std::cout << du_ << std::endl;
    // std::cout << du_*time_().deltaT << std::endl;
    // std::cout << u_ << std::endl;

    u_ += (du_ * time_.get_deltaT());

    log().info_end();
}

REGISTER_DEF_TYPE(TRANSPORTEQN, Momentum);
